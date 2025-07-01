function result = uwb_BasebandImpulseResponse(~, cfg, plotFlag)
% uwb_BasebandImpulseResponse - UWB基带冲激响应分析（符合HRP UWB物理层标准）
% 这是一个内部示例辅助函数，随时可能更改。
%
% 功能描述：
%   执行HRP UWB物理层发射机脉冲互相关函数检查，确保：
%   1. 主瓣幅值大于等于0.8
%   2. 主瓣幅度>=0.8的持续时间满足表16-28定义的要求（根据信道动态设置Tw）
%   3. 所有旁瓣幅值不超过0.3
%   根据信道编号自动设置对应的主瓣宽度Tw参数（参考IEEE 802.15.4a表16-28）
%
% 输入:
%   ~ - 波形数据（暂不使用，预留接口）
%   cfg - 配置参数结构体，必须包含：
%         .Channel - 信道编号（0,3,5,6,7,8,10,11,12,14,15等）
%         .CrossCorrMainPeakThreshold - 主峰幅值阈值（默认0.8）
%         .CrossCorrSideLobeThreshold - 旁瓣幅值阈值（默认0.3）
%         注意：MainLobeDurationMin会根据信道自动设置为对应的Tw值
%   plotFlag - 可选参数，是否绘制图形 (默认为false)
%
% 输出:
%   result - 结构体，包含检查结果:
%            .crossCorrPassed: 互相关检查是否通过 (logical)
%            .timeMaskPassed: 时域模板检查是否通过 (logical)
%            .overallPassed: 整体检查是否通过 (logical)
%            .crossCorrAnalysis: 详细的互相关分析结果
%            .message: 检查结果描述

if nargin < 3
    plotFlag = false; % 默认不绘图
end

%% 根据信道编号设置主瓣宽度Tw（参考表16-28）
% 根据IEEE 802.15.4a标准表16-28设置各信道对应的主瓣宽度参数
channelNumber = cfg.Channel;
switch channelNumber
    case {0, 3, 5, 6, 8, 10, 12, 14}  % 信道组{0;3;5;6;8;10;12;14}
        Tw_ns = 0.5;  % 主瓣宽度，纳秒
        Tp_ns = 2.00; % 脉冲持续时间，纳秒
    case 7
        Tw_ns = 0.2;  % 主瓣宽度，纳秒
        Tp_ns = 0.92; % 脉冲持续时间，纳秒
    case {4, 11}  % 信道组{4;11}
        Tw_ns = 0.2;  % 主瓣宽度，纳秒
        Tp_ns = 0.75; % 脉冲持续时间，纳秒
    case 15
        Tw_ns = 0.2;  % 主瓣宽度，纳秒
        Tp_ns = 0.74; % 脉冲持续时间，纳秒
    otherwise
        % 默认值（对于未定义的信道）
        Tw_ns = 0.5;  % 主瓣宽度，纳秒
        Tp_ns = 2.00; % 脉冲持续时间，纳秒
        warning('未识别的信道编号 %d，使用默认参数 Tw=%.1f ns, Tp=%.2f ns', ...
               channelNumber, Tw_ns, Tp_ns);
end

% 更新配置参数中的主瓣持续时间最小值
cfg.MainLobeDurationMin = Tw_ns;  % 使用表16-28中的Tw值作为最小主瓣宽度要求

fprintf('信道 %d: 脉冲持续时间 Tp=%.2f ns, 主瓣宽度 Tw=%.1f ns\n', ...
        channelNumber, Tp_ns, Tw_ns);

spc = cfg.SamplesPerPulse; % 每个脉冲的采样点数

%% 创建巴特沃斯脉冲

% 创建一个4阶巴特沃斯滤波器，3db带宽（截止频率）为500 MHz
N = cfg.ButterworthOrder; % 滤波器阶数
Fc = cfg.ButterworthCutoff; % 截止频率
Fs = Fc*spc; % 采样频率
[b,a] = uwb_butter(N, Fc/Fs); % 设计巴特沃斯滤波器

% 输入一个冲激信号到巴特沃斯滤波器，生成巴特沃斯脉冲
Tp = Tp_ns; % ns，脉冲宽度，根据表16-28和信道编号动态设置
NTp = cfg.PulseWidthMultiplier; % 脉冲宽度的倍数
xSpan = NTp*Tp; % ns，总时间跨度
impulse = [1; zeros((NTp*spc-1), 1)]; % 冲激信号
pulse = filter(b, a, impulse); % 通过滤波器得到脉冲

% 补偿滤波器延迟，使波形在t=0处居中
[~, idx] = max(pulse); % 找到脉冲最大值的位置
% butterPulse = circshift(butterPulse, round(length(butterPulse)/2)-idx);
len = length(pulse); % 脉冲长度
pulseCentered = [zeros(round(len/2)-idx, 1); pulse(1:(end-round(len/2)+idx))]; % 居中脉冲

% 互相关性检查:直接进行数据计算和检查

%% 创建根升余弦脉冲
beta = cfg.RRCBeta; % 滚降系数
t = -xSpan/2:(xSpan/(len-1)):xSpan/2; % 时间向量
tTp = t/Tp; % 归一化时间
r = (4*beta/pi*sqrt(Tp)) * (cos((1+beta)*pi*tTp) + sin((1-beta)*pi*tTp)./(4*beta*tTp))./(1-(4*beta*tTp).^2); % 根升余弦脉冲
r = r'; % 转置为列向量

% 避免NaN值:
r(t==0) = (1+beta*(4/pi -1))/Tp; % t=0处的特殊处理

% 归一化:
r = r * 1/max(r); % 最大值归一化

% 计算巴特沃斯脉冲与根升余弦脉冲的互相关
x = uwb_xcorr(r, pulseCentered, 'normalized'); % 归一化互相关

%% 执行互相关检查
% 互相关检查阈值（从配置参数读取，符合UWB标准）
T1 = cfg.CrossCorrMainPeakThreshold;  % 主峰检查阈值
T2 = cfg.CrossCorrSideLobeThreshold;  % 旁瓣检查阈值


% 检查主峰
crossCorrMainPeakOK = any(x(round(end/2 +[-1 0 1]))>T1); % 检查主峰,0延迟处任意三个点中一个大于0.8

% 旁瓣检查：检查所有大于T2的区域是否连续
ind = find(abs(x)>T2); % 找到所有幅值大于T2的点
if ~isempty(ind)
    % 检查这些点是否构成连续区域
    crossCorrSideLobesOK = isequal(ind', min(ind):max(ind));
    % 如果不连续，说明存在分离的旁瓣超过阈值
    if ~crossCorrSideLobesOK
        sideLobeViolationMsg = 'Cross-correlation greater than 0.3 is non-contiguous. That means some side-lobe is greater than 0.3';
    else
        sideLobeViolationMsg = '';
    end
else
    % 没有任何点超过T2，旁瓣检查通过
    crossCorrSideLobesOK = true;
    sideLobeViolationMsg = '';
end

% 主瓣幅度大于0.8的持续时间检查（符合HRP UWB物理层标准要求）
% 根据表16-28的定义，使用信道对应的主瓣宽度Tw作为最小持续时间要求
MIN_MAINLOBE_DURATION_NS = cfg.MainLobeDurationMin;  % 已根据信道设置为对应的Tw值 (ns)

% 找到主瓣区域（幅值>=0.8的连续区域）
mainLobeIndices = find(x >= T1);
if ~isempty(mainLobeIndices)
    % 确保主瓣区域是连续的
    diffIndices = diff(mainLobeIndices);
    if all(diffIndices == 1) || length(mainLobeIndices) == 1
        % 主瓣是连续的，计算幅度>=0.8的持续时间
        mainLobeDuration_samples = length(mainLobeIndices);
        % 计算时间分辨率：总时间跨度/(采样点数-1)
        timeResolution = (2*xSpan) / (length(x) - 1);  % ns per sample for correlation
        mainLobeDuration_ns = mainLobeDuration_samples * timeResolution;
        
        % 检查主瓣幅度>=0.8的持续时间是否符合要求（只检查最小值）
        mainLobeDurationOK = mainLobeDuration_ns >= MIN_MAINLOBE_DURATION_NS;
    else
        % 主瓣不连续，视为不符合要求
        mainLobeDuration_ns = 0;
        mainLobeDurationOK = false;
    end
else
    % 没有找到符合幅值要求的主瓣
    mainLobeDuration_ns = 0;
    mainLobeDurationOK = false;
end

% 综合互相关检查结果（包含主瓣幅度>=0.8的持续时间要求）
crossCorrPassed = crossCorrMainPeakOK && crossCorrSideLobesOK && mainLobeDurationOK;

%% 时域模板检查（完整的UWB脉冲时域特性分析）
% 脉冲归一化
pulseNormalized = pulse/max(abs(pulse));

% 1. 脉冲宽度测量（半功率点 -3dB）
pulseAmplitude = abs(pulseNormalized);
halfPowerLevel = max(pulseAmplitude) / sqrt(2); % -3dB点
halfPowerIndices = find(pulseAmplitude >= halfPowerLevel);
if ~isempty(halfPowerIndices)
    pulseWidth_ns = (halfPowerIndices(end) - halfPowerIndices(1) + 1) * (xSpan / (len-1));
else
    pulseWidth_ns = 0;
end

% 2. 上升时间测量（10%-90%幅度）
[maxAmp, maxIdx] = max(pulseAmplitude);
level10 = 0.1 * maxAmp;
level90 = 0.9 * maxAmp;

% 寻找上升沿的10%和90%点
risingEdge = pulseAmplitude(1:maxIdx);
idx10_rising = find(risingEdge >= level10, 1, 'first');
idx90_rising = find(risingEdge >= level90, 1, 'first');

if ~isempty(idx10_rising) && ~isempty(idx90_rising)
    riseTime_samples = idx90_rising - idx10_rising;
    riseTime_ns = riseTime_samples * (xSpan / (len-1));
else
    riseTime_ns = 0;
end

% 3. 下降时间测量（90%-10%幅度）
fallingEdge = pulseAmplitude(maxIdx:end);
idx90_falling = find(fallingEdge <= level90, 1, 'first') + maxIdx - 1;
idx10_falling = find(fallingEdge <= level10, 1, 'first') + maxIdx - 1;

if ~isempty(idx90_falling) && ~isempty(idx10_falling)
    fallTime_samples = idx10_falling - idx90_falling;
    fallTime_ns = fallTime_samples * (xSpan / (len-1));
else
    fallTime_ns = 0;
end

% 4. 过冲/欠冲检查
% 寻找脉冲前后的稳态值
preStableRegion = pulseNormalized(1:round(len*0.1));   % 前10%区域
postStableRegion = pulseNormalized(round(len*0.9):end); % 后10%区域
baselineLevel = mean([preStableRegion; postStableRegion]);

% 计算过冲（正方向超调）
overshoot = max(pulseNormalized) - 1; % 相对于归一化峰值的过冲
if overshoot < 0, overshoot = 0; end

% 计算欠冲（负方向最大偏移）
undershoot = min(pulseNormalized) - baselineLevel;
if undershoot > 0, undershoot = 0; end
undershoot = abs(undershoot);

% 5. 脉冲对称性检查
leftHalf = pulseAmplitude(1:maxIdx);
rightHalf = pulseAmplitude(maxIdx:end);
% 调整长度使其相等
minLen = min(length(leftHalf), length(rightHalf));
leftHalf = leftHalf(end-minLen+1:end);
rightHalf = rightHalf(1:minLen);
% 计算对称性相关系数
symmetryCorr = corrcoef(leftHalf, flipud(rightHalf));
if size(symmetryCorr, 1) > 1
    symmetryCoeff = symmetryCorr(1,2);
else
    symmetryCoeff = 1; % 完全对称
end

% 6. 频谱纯度检查（时域）
% 计算脉冲的THD（总谐波失真）近似
pulseFFT = fft(pulseNormalized, 2048);
pulsePSD = abs(pulseFFT).^2;
fundamentalIdx = find(pulsePSD == max(pulsePSD), 1);
harmonicPower = sum(pulsePSD) - pulsePSD(fundamentalIdx);
THD_approx = sqrt(harmonicPower / pulsePSD(fundamentalIdx)) * 100; % 百分比

%% 时域模板合规性检查
% UWB脉冲时域标准要求（IEEE 802.15.4a）
PULSE_WIDTH_MIN = 0.5;    % ns，最小脉冲宽度
PULSE_WIDTH_MAX = 10.0;   % ns，最大脉冲宽度
RISE_TIME_MAX = 2.0;      % ns，最大上升时间
FALL_TIME_MAX = 2.0;      % ns，最大下降时间
OVERSHOOT_MAX = 0.1;      % 10%，最大过冲
UNDERSHOOT_MAX = 0.1;     % 10%，最大欠冲
SYMMETRY_MIN = 0.8;       % 最小对称性系数
THD_MAX = 20;             % %，最大总谐波失真

% 执行各项检查
pulseWidthOK = (pulseWidth_ns >= PULSE_WIDTH_MIN) && (pulseWidth_ns <= PULSE_WIDTH_MAX);
riseTimeOK = riseTime_ns <= RISE_TIME_MAX;
fallTimeOK = fallTime_ns <= FALL_TIME_MAX;
overshootOK = overshoot <= OVERSHOOT_MAX;
undershootOK = undershoot <= UNDERSHOOT_MAX;
symmetryOK = symmetryCoeff >= SYMMETRY_MIN;
thdOK = THD_approx <= THD_MAX;

% 综合时域模板检查结果
timeMaskPassed = pulseWidthOK && riseTimeOK && fallTimeOK && ...
                 overshootOK && undershootOK && symmetryOK && thdOK;

% 准备返回结果
result.crossCorrPassed = crossCorrPassed;
result.timeMaskPassed = timeMaskPassed;
result.overallPassed = crossCorrPassed && timeMaskPassed;

% 添加互相关详细检查结果
result.crossCorrAnalysis.mainPeakOK = crossCorrMainPeakOK;
result.crossCorrAnalysis.sideLobesOK = crossCorrSideLobesOK;
result.crossCorrAnalysis.mainLobeDurationOK = mainLobeDurationOK;
result.crossCorrAnalysis.mainLobeDuration_ns = mainLobeDuration_ns;
result.crossCorrAnalysis.mainLobeThreshold = T1;
result.crossCorrAnalysis.sideLobeThreshold = T2;
result.crossCorrAnalysis.sideLobeViolationMsg = sideLobeViolationMsg;
% 添加Tw和信道信息
result.crossCorrAnalysis.mainLobeWidthTw_ns = Tw_ns;
result.crossCorrAnalysis.channelNumber = channelNumber;
result.crossCorrAnalysis.pulseDurationTp_ns = Tp_ns;

% 添加详细的时域分析结果
result.timeDomainAnalysis.pulseWidth_ns = pulseWidth_ns;
result.timeDomainAnalysis.riseTime_ns = riseTime_ns;
result.timeDomainAnalysis.fallTime_ns = fallTime_ns;
result.timeDomainAnalysis.overshoot = overshoot;
result.timeDomainAnalysis.undershoot = undershoot;
result.timeDomainAnalysis.symmetryCoeff = symmetryCoeff;
result.timeDomainAnalysis.THD_percent = THD_approx;

% 添加各项检查结果
result.timeDomainChecks.pulseWidthOK = pulseWidthOK;
result.timeDomainChecks.riseTimeOK = riseTimeOK;
result.timeDomainChecks.fallTimeOK = fallTimeOK;
result.timeDomainChecks.overshootOK = overshootOK;
result.timeDomainChecks.undershootOK = undershootOK;
result.timeDomainChecks.symmetryOK = symmetryOK;
result.timeDomainChecks.thdOK = thdOK;

% 输出检查结果
fprintf('=== UWB基带冲激响应分析结果 ===\n');
if crossCorrPassed
    fprintf('✓ 互相关检查通过\n');
    fprintf('  - 主峰检查: %s\n', char(crossCorrMainPeakOK*'通过' + ~crossCorrMainPeakOK*'失败'));
    fprintf('  - 旁瓣检查: %s (阈值: %.1f)\n', ...
            char(crossCorrSideLobesOK*'通过' + ~crossCorrSideLobesOK*'失败'), T2);
    fprintf('  - 主瓣幅度>=0.8持续时间检查: %s (%.3f ns, 要求Tw: ≥%.1f ns (信道%d))\n', ...
            char(mainLobeDurationOK*'通过' + ~mainLobeDurationOK*'失败'), ...
            mainLobeDuration_ns, MIN_MAINLOBE_DURATION_NS, channelNumber);
else
    fprintf('✗ 互相关检查失败\n');
    fprintf('  - 主峰检查: %s\n', char(crossCorrMainPeakOK*'通过' + ~crossCorrMainPeakOK*'失败'));
    fprintf('  - 旁瓣检查: %s (阈值: %.1f)\n', ...
            char(crossCorrSideLobesOK*'通过' + ~crossCorrSideLobesOK*'失败'), T2);
    if ~isempty(sideLobeViolationMsg)
        fprintf('    警告: %s\n', sideLobeViolationMsg);
    end
    fprintf('  - 主瓣幅度>=0.8持续时间检查: %s (%.3f ns, 要求Tw: ≥%.1f ns (信道%d))\n', ...
            char(mainLobeDurationOK*'通过' + ~mainLobeDurationOK*'失败'), ...
            mainLobeDuration_ns, MIN_MAINLOBE_DURATION_NS, channelNumber);
end

% 详细的时域模板检查输出
if timeMaskPassed
    fprintf('✓ 时域模板检查通过\n');
else
    fprintf('✗ 时域模板检查失败\n');
end

fprintf('  时域特性分析:\n');
fprintf('  - 脉冲宽度: %.3f ns (要求: %.1f-%.1f ns) %s\n', ...
        pulseWidth_ns, PULSE_WIDTH_MIN, PULSE_WIDTH_MAX, ...
        char(pulseWidthOK*10004 + ~pulseWidthOK*10060));
fprintf('  - 上升时间: %.3f ns (要求: ≤%.1f ns) %s\n', ...
        riseTime_ns, RISE_TIME_MAX, ...
        char(riseTimeOK*10004 + ~riseTimeOK*10060));
fprintf('  - 下降时间: %.3f ns (要求: ≤%.1f ns) %s\n', ...
        fallTime_ns, FALL_TIME_MAX, ...
        char(fallTimeOK*10004 + ~fallTimeOK*10060));
fprintf('  - 过冲: %.2f%% (要求: ≤%.0f%%) %s\n', ...
        overshoot*100, OVERSHOOT_MAX*100, ...
        char(overshootOK*10004 + ~overshootOK*10060));
fprintf('  - 欠冲: %.2f%% (要求: ≤%.0f%%) %s\n', ...
        undershoot*100, UNDERSHOOT_MAX*100, ...
        char(undershootOK*10004 + ~undershootOK*10060));
fprintf('  - 对称性: %.3f (要求: ≥%.1f) %s\n', ...
        symmetryCoeff, SYMMETRY_MIN, ...
        char(symmetryOK*10004 + ~symmetryOK*10060));
fprintf('  - THD: %.1f%% (要求: ≤%.0f%%) %s\n', ...
        THD_approx, THD_MAX, ...
        char(thdOK*10004 + ~thdOK*10060));

if result.overallPassed
    result.message = sprintf('所有HRP UWB物理层测试项均通过，脉冲符合标准要求 (主瓣幅度>=0.8持续时间: %.2f ns)', mainLobeDuration_ns);
    fprintf('✓ 整体分析: 所有HRP UWB物理层测试项均通过\n');
    fprintf('  符合HRP UWB发射机互相关要求：主瓣≥%.1f，旁瓣≤%.1f，持续时间Tw≥%.1f ns (信道%d)\n', ...
            T1, T2, MIN_MAINLOBE_DURATION_NS, channelNumber);
else
    result.message = '部分HRP UWB物理层测试项失败，请检查脉冲参数以符合标准要求';
    fprintf('✗ 整体分析: 部分HRP UWB物理层测试项失败\n');
end
fprintf('===============================\n');

%% 可选绘图
if plotFlag
    %% 互相关性检查图:
    figXCorr = figure; % 新建图形窗口
    
    title('互相关性一致性检查') % 设置标题
    
    subplot(1, 3, 1) % 创建1行3列的第1个子图
    plot(-xSpan/2:(xSpan/(len-1)):xSpan/2, pulseCentered); % 绘制居中后的巴特沃斯脉冲
    title('巴特沃斯脉冲') % 设置子图标题
    xlabel('时间 (ns)') % 设置x轴标签
    axis([-NTp NTp min(pulseCentered) max(pulseCentered)]) % 设置坐标轴范围

    subplot(1, 3, 2) % 创建第2个子图
    plot(t, r); % 绘制根升余弦脉冲
    title('RRC脉冲') % 设置子图标题
    axis([-NTp NTp min(r) max(r)]) % 设置坐标轴范围
    xlabel('时间 (ns)') % 设置x轴标签
    
    subplot(1, 3, 3) % 创建第3个子图
    len_corr = length(x); % 互相关长度
    plot(-xSpan:(2*xSpan/(len_corr-1)):xSpan, x); % 绘制互相关曲线
    
    hold on
    plot([-xSpan/2 xSpan/2], [T1 T1], 'r-', 'LineWidth', 2) % 绘制T1阈值线
    plot([-xSpan/2 xSpan/2], [T2 T2], 'r:', 'LineWidth', 1.5) % 绘制T2阈值线
    plot([-xSpan/2 xSpan/2], [-T2 -T2], 'r:', 'LineWidth', 1.5) % 绘制-T2阈值线
    
    % 高亮显示主瓣区域
    if ~isempty(mainLobeIndices)
        timeAxis = -xSpan:(2*xSpan/(len_corr-1)):xSpan;
        mainLobeTimeStart = timeAxis(mainLobeIndices(1));
        mainLobeTimeEnd = timeAxis(mainLobeIndices(end));
        fill([mainLobeTimeStart mainLobeTimeEnd mainLobeTimeEnd mainLobeTimeStart], ...
             [T1 T1 max(x)*1.1 max(x)*1.1], 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        text(0, max(x)*0.8, sprintf('主瓣幅度>=0.8持续时间: %.2f ns (要求Tw≥%.1f ns, 信道%d)', ...
                                    mainLobeDuration_ns, MIN_MAINLOBE_DURATION_NS, channelNumber), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'white');
    end
    
    title('互相关') % 设置子图标题
    axis([-NTp NTp -0.5 max(x)]) % 设置坐标轴范围
    xlabel('时间 (ns)') % 设置x轴标签
    
    % 将图例放在左下角，避免遮挡主瓣区域的波形
    lgd = legend('互相关', '主瓣阈值 (0.8)', '旁瓣阈值 (±0.3)', '', '主瓣幅度>0.8区域', ...
                 'Location', 'southwest', 'FontSize', 8);
    % 设置图例背景半透明，减少对波形的影响
    lgd.Color = [1 1 1 0.9]; % 白色背景，90%不透明度
    lgd.EdgeColor = [0.7 0.7 0.7]; % 浅灰色边框
    lgd.Box = 'on'; % 显示边框
    
    % 在图上显示检查结果（居中置顶）
    if crossCorrPassed
        sgtitle(sprintf('互相关检查: 通过 ✓ (主瓣幅度>=0.8持续时间: %.2f ns, 要求Tw≥%.1f ns, 信道%d)', ...
                       mainLobeDuration_ns, MIN_MAINLOBE_DURATION_NS, channelNumber), ...
                'Color', 'green', 'FontSize', 14, 'FontWeight', 'bold');
    else
        failureReason = '';
        if ~crossCorrMainPeakOK, failureReason = [failureReason '主峰 ']; end
        if ~crossCorrSideLobesOK, failureReason = [failureReason '旁瓣 ']; end
        if ~mainLobeDurationOK, failureReason = [failureReason '持续时间 ']; end
        sgtitle(sprintf('互相关检查: 失败 ✗ (失败项: %s, 要求Tw≥%.1f ns, 信道%d)', ...
                       failureReason, MIN_MAINLOBE_DURATION_NS, channelNumber), ...
                'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold');
    end

%% 时域模板:

% 绘制脉冲及其时域模板:
figMask = figure; % 新建图形窗口
pulse = pulse/max(pulse); % 脉冲归一化

% 时域模板固定显示参数（符合UWB标准）
TIME_MASK_X_MIN = -3;      % x轴最小值 (Tp单位)
TIME_MASK_X_MAX = 9;       % x轴最大值 (Tp单位)
TIME_MASK_Y_MIN = -0.8;    % y轴最小值
TIME_MASK_Y_MAX = 1.1;     % y轴最大值
TIME_MASK_ALPHA = 0.75;    % 模板透明度

xPulseStart = cfg.PulseStartTime; % Tp为单位的起始时间
xPulseEnd = cfg.PulseEndTime; % Tp为单位的结束时间
plot((xPulseStart:(xPulseEnd-xPulseStart)/(length(pulse)-1):xPulseEnd), pulse, 'b-o'); % 绘制脉冲
xMin = TIME_MASK_X_MIN; % x轴最小值
xMax = TIME_MASK_X_MAX; % x轴最大值
yMin = TIME_MASK_Y_MIN; % y轴最小值
yMax = TIME_MASK_Y_MAX; % y轴最大值
axis([xMin xMax yMin yMax]) % 设置坐标轴范围
t0 = 0; % t=0

% 绘制模板区域
a = TIME_MASK_ALPHA; % 透明度
darkRed = [200, 0, 0]/255; % 深红色
patch([xMin xMin t0-1.25 t0-1.25 t0+1 t0+1 xMax xMax], [yMax 0.015 0.015 1 1 0.3 0.3 yMax], darkRed, 'FaceAlpha', a); % 上模板
patch([xMin xMin t0 t0 t0+2 t0+2 xMax xMax], [yMin -0.015 -0.015 -0.5 -0.5 -0.3 -0.3 yMin], darkRed, 'FaceAlpha', a); % 下模板

xlabel('时间 (Tp)') % 设置x轴标签
title('巴特沃斯脉冲的时域模板') % 设置标题

% 在图上显示检查结果（放在图形右上角）
if timeMaskPassed
    text(xMax*0.7, yMax*0.8, '时域模板检查: 通过 ✓', 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
else
    text(xMax*0.7, yMax*0.8, '时域模板检查: 失败 ✗', 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
end

assignin('caller', 'figXCorrelation', figXCorr); % 将互相关图句柄赋值到工作区
assignin('caller', 'figTimeDomainMask', figMask); % 将时域模板图句柄赋值到工作区

end

end

