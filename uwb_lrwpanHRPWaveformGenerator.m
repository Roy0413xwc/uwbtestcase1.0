% 是实现 IEEE 802.15.4a/z 标准 HRP（高脉冲重复频率）超宽带波形生成的核心函数。
% 该函数根据配置参数生成符合标准的 UWB 波形，支持 HPRF、BPRF 和传统 802.15.4a 三种模式。

function [wave, symbols] = uwb_lrwpanHRPWaveformGenerator(PSDU, cfg)
% 输入参数：
% PSDU：物理层协议数据单元（来自 MAC 层的二进制数据）
% cfg：配置对象（lrwpanHRPConfig类实例，包含模式、数据速率等参数）
% 输出参数：
% wave：生成的 UWB 波形（经过脉冲成形的模拟信号）
% symbols：调制后的符号序列（未进行脉冲成形的数字信号）

%LRWPANHRPWAVEFORMGENERATOR 生成符合IEEE 802.15.4a/z标准的HRP超宽带波形
%   WAVE = LRWPANHRPWAVEFORMGENERATOR(PSDU, CFG) 按照IEEE 802.15.4a/z标准创建HRP(高脉冲重复频率)超宽带波形WAVE。
%   CFG必须是一个<a href="matlab:help('lrwpanHRPConfig')">lrwpanHRPConfig</a>对象。
%   根据输入配置对象CFG的Mode属性，WAVE可以是HRPF(更高脉冲重复频率)、BPRF(基本脉冲重复频率)或IEEE 802.15.4a波形。
%
%   PSDU是物理层协议数据单元，由上层(MAC层)提供。
%   其长度必须是8的整数倍，且对于BPRF和802.15.4a模式，长度不能超过1016；对于HPRF模式，长度不能超过32760。
%   允许空PSDU([])。在所有情况下，PSDU的长度必须与CFG配置对象中的PSDULength属性值匹配。
%
%   [WAVE, SYMBOLS] = LRWPANHRPWAVEFORMGENERATOR(PSDU, CFG) 还会返回脉冲整形前的HRP信号，
%   即在调制(符号映射)和前导码插入之后的信号。

%#codegen

uwb_validate_Params(cfg); %调用配置验证函数，用于确保输入的配置对象合法

%% 验证PSDU长度
if strcmp(cfg.Mode, 'HPRF')   %为HPRF,相同返回1(true)，不同返回0
  numLenBits = 12;
else % BPRF, 15.4a
  numLenBits = 7;
end

PSDUlen = length(PSDU);
if PSDUlen > 8*(2^numLenBits-1) % 超过最大长度，报错
  error(message('lrwpan:LRWPAN:InvalidMPDULength', 8*(2^numLenBits-1)));
end
if PSDUlen ~= 8*cfg.PSDULength % 修改PSDUlen或PSDU，确保载荷长度和PSDUlen相同
  coder.internal.warning('lrwpan:LRWPAN:PSDULenMismatchTX', PSDUlen, 8*cfg.PSDULength);
  PSDUlen = min(PSDUlen, 8*cfg.PSDULength);
  PSDU = PSDU(1:PSDUlen);
end

%% 编码
% 在802.15.4a模式或STS包配置不为3时，需要传输有效载荷(payload)和物理层头部(PHR)
% In BPRF/HPRF modes, payload and PHR can be omitted and only STS transmitted
% 在 BPRF（基本脉冲重复频率）/HPRF（高脉冲重复频率）模式下，可省略有效载荷（Payload）和 PHR（物理层头部），仅传输 STS（超帧起始序列）
if strcmp(cfg.Mode, '802.15.4a') || cfg.STSPacketConfiguration ~= 3
  %% 1. Reed-Solomon 编码
  if ~strcmp(cfg.Mode, 'HPRF') || cfg.ConstraintLength ~= 7 
    %不为HPRF或卷积编码长度不为7进行RS编码
    encoding = true;
    rsPSDU = uwb_hrpRS(PSDU, encoding);
  else
    % 在HPRF模式下且卷积约束长度为7时不启用RS编码，直接输出
    rsPSDU = PSDU;
  end

  %% 2. PHY Header and SECDED encoding
  % 生成物理层帧头并对帧头进行单纠错双检测编码
  PHR = uwb_createPHRwithSECDED(PSDUlen, cfg);

  %% 3. Convolutional Encoding
  % 卷积编码
  convolCW = uwb_convolEnc(PHR, rsPSDU, cfg);

  %% 4. Symbol Mapper (Modulation)
  % 符号映射
  symbols = uwb_symbolMapper(convolCW, cfg);
else
  % 在 BPRF/HPRF模式下，STS包配置=3时，可省略有效载荷（Payload）和 PHR（物理层头部），仅传输 STS（超帧起始序列）
  % 传输符号为空
  symbols = [];
end

%% 5. Preamble Insertion
% 插入前导码，前导码由SYNC和SFD组成
SHR = uwb_createSHR(cfg);

%% 6. STS (Srambled Timestamp Sequence)
if ~strcmp(cfg.Mode, '802.15.4a') && cfg.STSPacketConfiguration ~= 0
  % STS 仅适用于 BPRF/HPRF 模式；当 STS 配置为 1、2 或 3 时
  % STS配置为1，STS序列长度最短适合低功耗场景，减少传输时间和功耗
  % STS配置为2，用于常规通信，平衡同步效率和硬件复杂度
  % STS配置为3，用于高可靠性需求，尤其在多径干扰严重的环境中
  
  % 生成STS序列
  STS = uwb_createSTS(cfg);  
  
  %根据不同的配置STS插入不同位置
  switch cfg.STSPacketConfiguration
    case 1
      symbols = [SHR; STS; symbols];
    case 2
      symbols = [SHR; symbols; STS];
    case 3
      symbols = [SHR; STS];
  end
else
  % No STS
  symbols = [SHR; symbols];
end

% 生成波形，经过butter滤波器，并转化为复数形式
wave = complex(uwb_butterworthFilter(symbols, cfg.SamplesPerPulse));


%% Subfunctions:部分fountion函数

%% 2. PHY Header and SECDED encoding
function phr = uwb_createPHRwithSECDED(PSDUlen, cfg)
% 参考2024版标准16.2.7  
% As per Sec. 15.2.7 in IEEE Std 802.15.4™‐2020

HPRF = strcmp(cfg.Mode, 'HPRF'); %相同的话，返回值为1

phr = zeros(13, 1);

if ~HPRF %若不为HPRF模式

  % 1. 0-1前两位为Data Rate，根据DataRateNum和MeanPRFNum设置
  rate = cfg.DataRateNum;
  if rate == 110
    phr(1:2) = [0; 0];
  elseif rate == 850
    phr(1:2) = [0; 1];
  elseif (rate == 6810 && cfg.MeanPRFNum ~= 3.9) || (rate == 1700 && cfg.MeanPRFNum == 3.9)
    phr(1:2) = [1; 0];
  elseif (rate == 27240 && cfg.MeanPRFNum ~= 3.9) || (rate == 6810 && cfg.MeanPRFNum == 3.9)
    phr(1:2) = [1; 1];
  % 否则，由validateConfig抛出错误，其他条件下不存在
  end

  % 2. 2-8 Frame Length 占七位
  len = PSDUlen/8;                 %将bit数转化为字节数
  phr(3:9) = int2bit(len, 7);      %函数用于将整数转换为二进制位表示，超过七位高位截断，不足七位，左侧补零

  % 3. 9 Ranging 测距标识符
  phr(10) = double(cfg.Ranging);

  % 4. 10 Reserved  保留位
  phr(11) = 0;

  % 5. 11-12 Preamble Duration
  preambleLen = cfg.PreambleDuration; %扩频前导SYNC码的重复次数
  if preambleLen == 16
    phr(12:13) = [0; 0];
  elseif preambleLen == 64
    phr(12:13) = [0; 1];
  elseif preambleLen == 1024
    phr(12:13) = [1; 0];
  elseif preambleLen == 4096
    phr(12:13) = [1; 1];
  % 有且仅有以上四种情况，其他情况不允许
  end
  
else % HPRF  
  
  %0-1 前两位的配置
  len = PSDUlen/8;
  if cfg.STSPacketConfiguration==2 && (cfg.ExtraSTSGapIndex>0 || cfg.ExtraSTSGapLength>0)
    % 在STS配置=2，加入额外的STS间隔
    phr(1:2) = int2bit(cfg.ExtraSTSGapIndex, 2);
  else
    % 其他配置下，作为序列长度的延伸
    phr(1) = bitget(len, 12);  % A1  提取整数 len 的第 12 位二进制值（从右向左数，最低位为第 1 位）
    phr(2) = bitget(len, 11);  % A0
  end
  
  
  % 2. Payload length
  % 调用 bitget 函数获取 10 个最低有效位（LSB右边是高位）, 
  % fliplr for left-msb（左边是高位）
  phr(3:12) = fliplr(bitget(len, 1:10));
  
  % 3. Ranging:
  phr(13) = double(cfg.Ranging);
end

% END. Hamming coding - SECDED, i.e., Single error correction, double error detection
% 生成帧头之后对其进行单纠错双检测编码
phr = uwb_hrpSECDED(phr);

%% 3. Convolutional Encoding
function convolCW = uwb_convolEnc(PHR, rsPSDU, cfg)
% 参考标准16.3.3

% Two zeros for ConstraintLength = 3, six for length = 7
tailField = zeros(cfg.ConstraintLength-1, 1);
% 尾比特将卷积编码器重置为0，确保解码时状态同步

if cfg.ConvolutionalCoding   % 启用1/2维特比速率
  if ~strcmp(cfg.Mode, 'HPRF') || cfg.ConstraintLength == 3  % 不为HPRF或者卷积长度为3，
    % 具体原理见16.3.3.3
    % 见翻译标准表16-2
    convolIn = [PHR; rsPSDU; tailField];
    
  else % CL = 7
    % 见表16-2，添加尾比特，分别在 PHR（物理层头部）和 PSDU（物理层服务数据单元）后追加 6 个零比特
    convolIn = [PHR; tailField; rsPSDU; tailField];
  end
else  % 启用1维特比速率   不进行卷积编码
  % 见表16-2
  % 有一些情况不对PSDU进行编码，见翻译标准16-4
  % 一直对PHR进行卷积编码
  convolIn = [PHR; tailField];

end

% Rate 1/2 coding:速率为1/2
if ~(strcmp(cfg.Mode, 'HPRF') && cfg.ConstraintLength == 7)
  % 见16.3.3.3 
  trellis3 = poly2trellis(3, [2 5]);
  convolCW = convenc(convolIn, trellis3);
else
  % 见16.3.3.3
  trellis7 = poly2trellis(7, [133 171]);
  convolCW = convenc(convolIn, trellis7);  % repeat for codegen
end

if ~cfg.ConvolutionalCoding
  % Table 16-3 
  % 不启用卷积编码，仅对PHR编码，PSDU直接通过
  convolCW = [convolCW; rsPSDU];
end


%% 4. Symbol Mapper (Modulation)
function symbols = uwb_symbolMapper(convolCW, cfg)
% 非 HPRF 模式：基于跳时扩频（TH-PPM），通过校验位调制极性，系统位控制突发位置。
% HPRF 模式：使用查表映射和扩频加扰，适用于高脉冲重复频率场景。
% As per Sec. 15.3 in IEEE Std 802.15.4™‐2020

persistent pn % Persistent comm.PNSequence, as one-time setup is the computational bottleneck
% 声明持久变量，在函数多次调用之间保留变量值，避免重复初始化

% 0. Repackage input with 1 codeword per row
% 重新打包输入数据，使每行包含一个码字
tmp = reshape(convolCW, 2, []);% 两行
cws = tmp';% 两列 每行为一个码字（2比特：系统位+校验位）
numSym = size(cws, 1);% 获取符号总数，即两个比特为一个符号
inHPRF = strcmp(cfg.Mode, 'HPRF');

phrLen = 19;% 帧头原始长度
if ~inHPRF || cfg.ConstraintLength == 3
  numPHRsym = phrLen+2;
else % CL = 7 
  numPHRsym = phrLen+6;
end  %根据调制模式和卷积约束长度调整 PHR 符号数

% 2. Create spreading sequence obj
% 生成扩频序列
code = uwb_HRPCodes(cfg.CodeIndex);%获取扩频码
codeHat = code;
% a. Remove all zeros
codeHat(codeHat==0) = [];% 将0变为-1
% b. Replace all -1 with zeros
codeHat(codeHat==-1) = 0;
% c. Keep the 1st 15
codeHat = codeHat(1:15);% 取前15个元素
initialConditions = fliplr(codeHat); % 反转作为初始条件  % fliplr because of the difference in convention
% 生成 PN 序列的初始条件，用于后续扩频

% do a single PN sequence call, to facilitate codegen (pass InitialConditions and SamplesPerFrame as inputs)
% 计算每个符号的码片数
if ~inHPRF
  numChips = cfg.ChipsPerBurst;% 非 HPRF 模式：使用突发（Burst）结构，每个符号包含多个突发。
  % 这是一个突发里的chip数，需要两级配置
  % 只需要调制突发那一段，计算需要调制的chip数
else
  numChips = cfg.ChipsPerSymbol; % HPRF 模式：根据 PRF（脉冲重复频率）计算符号速率
  % chip 可能是标量或向量
  % PHR 和 PSDU 使用不同码片数（如ChipsPerBurst = [16, 8]）
end
totalSamples = numPHRsym*numChips(1) + (numSym-numPHRsym)*numChips(end);
% 总chip数 = PHR符号数 × PHR符号的码片数 + PSDU符号数 × PSDU符号的码片数
% The maximum value of totalSamples under all configurations is about 600K.
maxSamples = 6.1e5;


% 标准与 comm.PNSequence 的表示法存在差异。
% 标准中采用额外的连接（D^14），该连接靠近循环反馈的寄存器，
% 而 comm.PNSequence 中对应的连接为 D^1
if isempty(pn)
  pn = comm.PNSequence(Polynomial = '1 + D + D15', ...
    InitialConditionsSource = 'Input port', ...
    VariableSizeOutput=true, ...
    MaximumOutputSize=[maxSamples 1], ...
    Mask = 15); % Skip the first 15 (initial conditions), see example in Table 15-10
end % 使用 15 阶多项式生成 PN 序列，初始条件来自扩频码处理结果。
% Mask = 15表示跳过前 15 个样本（与标准中表 15-10 示例一致）

allPNSamples = pn(initialConditions, totalSamples);
reset(pn);  % prepare for next-function call

% 3. 调制每个码字 / 符号
% 3a. 前 21 个符号（PHR）的调制速率最高为 850 kb/s，即每个脉冲串至少 16 个码片
% 3b. 剩余符号（numSym-21）采用 cfgObj.DataRate 参数指定的速率调制
samplesPerBurst = nan; % % 为代码生成初始化变量，确保所有代码路径均有定义
if ~inHPRF
  % sps每个符号的chip数
  samplesPerBurst = cfg.ChipsPerBurst;
  sps = cfg.BurstsPerSymbol*samplesPerBurst;%每个符号的突发数乘每个突发的样本数就是整个样本
  % 非 HPRF 模式：使用突发（Burst）结构，每个符号包含多个突发。
else
  sps = cfg.ChipsPerSymbol*2*249.6/cfg.MeanPRFNum; % *2的目的感觉像加上保护间隔
  % 249.6是 IEEE 802.15.4z HPRF 模式的基准 PRF 单位，单位是 MHz
  % 表达式 2 * 249.6 / cfg.MeanPRFNum 是一个时间/速率归一化因子，
  % 用于将实际 PRF 下的采样长度，换算为标准 PRF（249.6 MHz）条件下的“等效采样点数”。
  % HPRF直接设置每个符号的扩频码片数
  % HPRF 模式：根据 PRF（脉冲重复频率）计算符号速率。
end

% sps是1或2元素向量。若为2元素，第1个值用于PHR，第2个值用于PSDU
PHRorPSDU = 1;
PHR_ID = 1;
PSDU_ID = length(sps);
symbols = zeros(numPHRsym*sps(PHR_ID)+(numSym-numPHRsym)*sps(PSDU_ID), 1);
% 是在初始化整个物理层（PHY）符号序列对应的采样点数组，也就是为整帧信号分配空间，单位是采样点数。

for symIdx = 1:numSym % 每个符号调制，主循环便利每一个符号
  % 处理PHR到PSDU的过渡（如果速率不同）
  % PHR（物理层头部）：通常采用较低调制速率（如 850 kb/s）以确保可靠性。
  % PSDU（有效载荷）：可能采用更高速率（如 1.6/6.8 Mb/s）以提升吞吐量。
  if ((inHPRF && cfg.ConstraintLength == 3) || (~inHPRF && ~isscalar(cfg.ChipsPerBurst))) && symIdx == numPHRsym+1
    PHRorPSDU = 1+double(symIdx>numPHRsym);       % isscalar用于判断一个变量是不是标量，是标量，返回1
  end % PHRorPSDU 被置成 2，把后续处理逻辑（采样率、扩频参数、编码参数、交织深度等）切换到 PSDU 的专用配置
  
  
  % 逐个获取每个码字，并逐个构建每个符号
  if PHRorPSDU == 1
    offset = 0;
    currSym = symIdx;
  else
    offset = numPHRsym*numChips(1);% 用于从 allPNSamples 中偏移跳过 PHR 段的扩频序列
    currSym = symIdx-numPHRsym;
  end

  % 区分 PHR（帧头）和 PSDU（有效载荷）部分，可能使用不同的调制参数
  % 获取当前符号的扩频序列
  spreadingSeq  = allPNSamples(offset+1+(currSym-1)*numChips(PHRorPSDU): offset+currSym*numChips(PHRorPSDU), 1);
  systematicBit = cws(symIdx, 1);% 系统位
  parityBit     = cws(symIdx, 2);% 极性位
  
  thisSymbol = zeros(sps(PHRorPSDU), 1);

  if ~inHPRF
    % 3. Calculate burst hopping position
    % 计算跳时位置
    m = log2(cfg.NumHopBursts);

    % 标准规定：当 m > Ncpb 时，LFSR 的时钟次数不得超过每突发码片数（ChipsPerBurst）
    m = min(m, cfg.ChipsPerBurst(PHRorPSDU));  % 限制跳变位数不超过每个突发的码片数
    burstPos = (2.^(0:m-1))*spreadingSeq(1:m);
    % 将比特转为整数索引，表示突发将放在 [0, NumHopBursts-1] 中的哪一位
    % burstPos ∈ [0, NumHopBursts - 1]，是一个可控范围

    % 4. Create burst
    % 生成突发信号
    thisBurst = (1-2*parityBit)*(1-2*spreadingSeq');
    thisBurst = thisBurst(:);
                                                                                                                                                                            
    % 5. Place burst in the corresponding position
    % 放置突发到指定位置
    burstNo = burstPos(1) + systematicBit*cfg.BurstsPerSymbol/2; % 使用 0 基索引；此外，burstPos 为标量，但添加 (1) 有助于代码生成
    thisSymbol(1+burstNo*samplesPerBurst(PHRorPSDU):(burstNo+1)*samplesPerBurst(PHRorPSDU)) = thisBurst;
    
  else % HPRF
    len = 4*249.6/cfg.MeanPRFNum;
    % IEEE 802.15.4z HPRF 模式中定义的 基本频率单元（MHz）
    symbolMap = uwb_hrpHPRFSymbolMap(cfg.MeanPRFNum, cfg.ConstraintLength);
    
    if PHRorPSDU == 1 % PHR
      numBits = size(symbolMap, 2);
    else % PSDU
      numBits = 2*len;% 是当前符号或数据单元中有效比特的数量，每个符号两个比特
    end % numBits调制的总比特数
    
    % 根据比特组合查表映射
    thisMapping = symbolMap(1+bit2int([systematicBit parityBit]', 2, false), 1:numBits);
    scrambled = (1-2*thisMapping) .* (1-2*spreadingSeq');
    
    % 加扰处理

    % 添加保护带和零填充
    scrambled = [reshape(scrambled, len, []); ...
                 zeros(len, numBits/len)];
    scrambled = scrambled(:);
    
    if cfg.MeanPRFNum == 124.8  
        % 相当于每一个原始比特之间插一个 0，进一步稀疏化符号，用于防止相邻脉冲干扰，也方便更精准的脉冲定位。
        scrambled = [scrambled'; 
                     zeros(1, numBits*2)];
        scrambled = scrambled(:);
    end
    thisSymbol = scrambled';
    thisSymbol = thisSymbol(:);
  end
  
  % 将每个符号放置到输出序列的正确位置，区分 PHR 和 PSDU 部分
  if PHRorPSDU == 1  % PHR
    startSymbolPos = 1+(symIdx-1)*sps(PHR_ID);
  else % PSDU
    phrEnd = numPHRsym*sps(PHR_ID);
    startSymbolPos = phrEnd + 1 +(symIdx-1-numPHRsym)*sps(PSDU_ID);
  end
  endSymbolPos = startSymbolPos + sps(PHRorPSDU)-1;
  symbols(startSymbolPos:endSymbolPos) = thisSymbol;
end



%% 5. Preamble Insertion
function SHR = uwb_createSHR(cfg)

%% SYNC
% As per Sec. 15.2.6.2 in IEEE Std 802.15.4™‐2020
code = uwb_HRPCodes(cfg.CodeIndex);
L = cfg.PreambleSpreadingFactor;  %前导码扩频因子
N = cfg.PreambleDuration;  %前导码重复次数，对应R&S软件中的sync length

% 1. 扩频处理：在每个三元符号(-1,0,+1)后插入L-1个0
spread = zeros(L*length(code), 1);
spread(1:L:end) = code;% 每隔L个位置插入一个码符号，其余为0

% 2. 生成同步序列(SYNC)：将扩频序列重复N次
SYNC = repmat(spread, N, 1); %列向量


%% SFD (Start of Frame delimiter)
% 3. 生成帧起始定界符(SFD)
% 见标准 16.2.6.3 
seq = uwb_getSFD(cfg);%获取SFD序列
SFD = seq.*spread; %SFD与扩频因子相乘
SFD = SFD(:);  %将SFD转为列向量
% 4. 组合SYNC和SFD形成SHR
SHR = [SYNC; SFD];


%% %% 6. STS (Srambled Timestamp Sequence)
function STS = uwb_createSTS(cfg)
% createSTS函数生成 STS 序列，用于 UWB 测距和定位中的时间戳标记，
% 支持 HPRF (高速率) 和 BPRF (基本速率) 两种模式
% 见标准16.2.9

len512 = 512; % chips   以512码片位单位  基本时间间隔
gap = zeros(len512, 1);
STS = gap; % 所有STS以间隔为开始

if strcmp(cfg.Mode, 'HPRF') && cfg.STSPacketConfiguration == 2
  % 在STS分组配置值=2的情况下，在PSDU与STS之间插入一个可选的附加间隔段
  STS = [STS; zeros(4*cfg.ExtraSTSGapLength, 1)];
end

if strcmp(cfg.Mode, 'BPRF') % BPRF
  numSegments = 1;
  segLen = 64*512; % STSSegmentLength的单位是512码片
  spreadingF = 8; % 扩频因子
else % HPRF
  numSegments = cfg.NumSTSSegments;
  segLen = cfg.STSSegmentLength*512; % STSSegmentLength的单位是512码片
  % 段长度由配置决定
  spreadingF = 4;
end

singleDRBGlen = 128; % 单个DRBG序列长度(比特)

% 获取已为以下参数预生成的 DRBG 随机比特：
% 密钥 ='4a5572bc90798c8e518d2449092f1b55'，
% 高 96 位 ='68debd3a599939dd57fdbb0e'，
% 低 32 位 ='2a10fac0'
% 这些值以 uint8 格式打包，需转换为比特：

% 加载预先生成的DRBG序列（用于加扰）
s = coder.load('allDRBG_STS.mat');
allDRBG = s.allDRBG;

counter = 1;

for idx = 1:numSegments
  activeSTS = zeros(segLen, 1);
  
  numDBRG = (segLen/(singleDRBGlen*spreadingF));% 每个段需要的DRBG块数
  for idx2 = 1:numDBRG
    % 读取DRBG序列并转换为比特
    singleDRBG_uint8 = allDRBG(counter, :);
    tmp = int2bit(singleDRBG_uint8, 8);
    singleDRBG = double(tmp(:));
    
    % 比特映射：0→+1，1→-1（符合UWB信号表示）
    singleDRBG(singleDRBG==1) = -1;
    singleDRBG(singleDRBG==0) = +1;
    
    % 扩频处理：将每个比特扩展为spreadingF个码片（4或8）   
    spreadBits = [singleDRBG'; zeros(spreadingF-1, singleDRBGlen)];
    spreadBits = spreadBits(:);

    % 填充到当前STS段
    activeSTS(1+(idx2-1)*singleDRBGlen*spreadingF : idx2*singleDRBGlen*spreadingF) = spreadBits;

    counter = counter+1;
  end
  STS = [STS; activeSTS; gap]; %#ok<AGROW>
end


%% 7. Pulse shaping
function wave = uwb_butterworthFilter(symbols, spc)

% 1. Create a 4th-order Butterworth filter with 3-db bandwidth (cutoff) at 500 MHz
N = 4;       % 阶数
Fc = 500e6;  % 截止频率
Fs = Fc*spc; % 采样率
[b,a] = butter(N, Fc/Fs);

% 2. Pass impulses to the Butterworth filter to create Butterworth pulses
% 将脉冲信号通过巴特沃斯滤波器，生成巴特沃斯脉冲
impulses = zeros(length(symbols)*spc, 1);%根据采样率将原始序列扩展
impulses(1:spc:end) = symbols;
wave = filter(b, a, impulses);
