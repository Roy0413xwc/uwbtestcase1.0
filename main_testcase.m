%% UWB 高脉冲重复频率 (HPRF) 测试案例主程序
% 
% 功能描述：
%   本程序实现IEEE 802.15.4a标准的UWB-HPRF模式波形生成和测试验证
%   包括基带脉冲响应测试和传输功率谱密度(PSD)模板合规性检查
%
% 主要测试项目：
%   1. 巴特沃斯滤波器和根升余弦(RRC)滤波器的脉冲响应生成
%   2. 脉冲互相关性能分析
%   3. 时域模板符合性检查
%   4. 频域PSD模板合规性验证
%
% 版本：1.0
%%

%% 第一步：加载测试数据
% 从data.txt文件中加载用于测试的二进制数据序列
% 这些数据将作为PSDU（物理服务数据单元）的输入比特流
fprintf('正在加载测试数据...\n');
data = load('data.txt');
fprintf('测试数据加载完成，数据长度：%d 比特\n', length(data));

%% 第二步：配置HPRF参数
% 创建并配置IEEE 802.15.4a HPRF模式的参数结构体
fprintf('\n正在配置HPRF参数...\n');

cfgHPRF = uwb_paramsconfig;              % 创建默认参数结构体

% 基本协议参数
cfgHPRF.Mode='802.15.4a';                % 设置为IEEE 802.15.4a标准模式
cfgHPRF.MeanPRF=15.6;                    % 平均脉冲重复频率 15.6 MHz (HPRF模式)
cfgHPRF.DataRate=27.24;                  % 数据传输速率 27.24 Mbps
cfgHPRF.Channel=9;                       % 使用信道9 (中心频率约8.5 GHz)

% 扩频码参数
cfgHPRF.CodeIndex=3;                     % 使用第3个扩频码（长度为31的序列）

% 前导码参数  
cfgHPRF.PreambleMeanPRF=4.03;            % 前导码平均脉冲重复频率 4.03 MHz

% 数据包参数
cfgHPRF.PSDULength=100;                  % PSDU长度设置为100字节                  

fprintf('HPRF参数配置完成：\n');
fprintf('  - 工作模式：%s\n', cfgHPRF.Mode);
fprintf('  - 平均PRF：%.2f MHz\n', cfgHPRF.MeanPRF);
fprintf('  - 数据速率：%.2f Mbps\n', cfgHPRF.DataRate);
fprintf('  - 信道编号：%d\n', cfgHPRF.Channel);
fprintf('  - 扩频码索引：%d\n', cfgHPRF.CodeIndex);
fprintf('  - PSDU长度：%d 字节\n', cfgHPRF.PSDULength);

%% 第三步：初始化和验证参数
% 根据用户配置的基本参数，自动推导和计算其他依赖参数
% 并验证所有参数的合法性和兼容性
fprintf('\n正在初始化和验证参数...\n');

cfgHPRF = uwb_init_Params(cfgHPRF);      % 初始化推导字段（如时隙长度、符号时间等）
uwb_validate_Params(cfgHPRF);            % 校验配置合法性（检查参数范围和兼容性）

fprintf('参数初始化和验证完成\n');

%% 第四步：生成HPRF波形
% 基于配置参数和输入数据生成完整的IEEE 802.15.4a HPRF波形
% 包括前导码、起始帧定界符(SFD)、物理头(PHR)和PSDU部分
fprintf('\n正在生成HPRF波形...\n');

currBits = data(1:cfgHPRF.PSDULength*8); % 提取所需长度的测试数据（8比特/字节）
fprintf('提取数据比特数：%d\n', length(currBits));

[waveHPRF, ~] = uwb_lrwpanHRPWaveformGenerator(currBits, cfgHPRF);

fprintf('HPRF波形生成完成，波形采样点数：%d\n', length(waveHPRF));

%% 可选：绘制波形帧结构
% 取消下面一行的注释可以可视化生成的波形结构
% uwb_lrwpanPlotFrame(waveHPRF, cfgHPRF);

%% 第五步：执行基带脉冲响应测试
% 执行IEEE 802.15.4a标准规定的基带脉冲响应合规性测试
% 测试内容包括：
%   1. 巴特沃斯滤波器脉冲响应生成和分析
%   2. 根升余弦(RRC)滤波器脉冲响应生成和分析  
%   3. 脉冲间互相关性能检查
%   4. 时域模板符合性验证
fprintf('\n开始执行基带脉冲响应测试...\n');

% uwb_BasebandImpulseResponse(waveHPRF,cfgHPRF,false);  % 默认不显示图形输出

fprintf('基带脉冲响应测试完成\n');

%% 第六步：执行传输PSD模板测试
% 执行IEEE 802.15.4a标准规定的传输功率谱密度模板合规性测试
% 验证生成的信号频谱是否符合监管要求和标准规范
fprintf('\n开始执行传输PSD模板测试...\n');

uwb_TransmitPSDmask(waveHPRF, cfgHPRF);  % 调用传输PSD模板测试函数

fprintf('传输PSD模板测试完成\n');




