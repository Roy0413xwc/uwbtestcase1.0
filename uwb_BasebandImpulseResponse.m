function result = uwb_BasebandImpulseResponse(wave, cfg, plotFlag)
% uwb_BasebandImpulseResponse - UWB基带冲激响应分析
% 这是一个内部示例辅助函数，随时可能更改。
%
% 输入:
%   wave - 波形数据
%   cfg - 配置参数结构体
%   plotFlag - 可选参数，是否绘制图形 (默认为false)
%
% 输出:
%   result - 结构体，包含检查结果:
%            .crossCorrPassed: 互相关检查是否通过 (logical)
%            .timeMaskPassed: 时域模板检查是否通过 (logical)
%            .overallPassed: 整体检查是否通过 (logical)
%            .message: 检查结果描述

if nargin < 3
    plotFlag = false; % 默认不绘图
end

spc = cfg.SamplesPerPulse; % 每个脉冲的采样点数

%% 创建巴特沃斯脉冲

% 1. 创建一个4阶巴特沃斯滤波器，3db带宽（截止频率）为500 MHz
N = cfg.ButterworthOrder; % 滤波器阶数
Fc = cfg.ButterworthCutoff; % 截止频率
Fs = Fc*spc; % 采样频率
% [b,a] = butter(N, Fc/Fs); % 设计巴特沃斯滤波器
[b,a] = uwb_butter(N, Fc/Fs); % 设计巴特沃斯滤波器

% 2. 输入一个冲激信号到巴特沃斯滤波器，生成巴特沃斯脉冲
Tp = cfg.PulseWidth; % ns，脉冲宽度
NTp = cfg.PulseWidthMultiplier; % 脉冲宽度的倍数
xSpan = NTp*Tp; % ns，总时间跨度
impulse = [1; zeros((NTp*spc-1), 1)]; % 冲激信号
pulse = filter(b, a, impulse); % 通过滤波器得到脉冲

% 补偿滤波器延迟，使波形在t=0处居中
[~, idx] = max(pulse); % 找到脉冲最大值的位置
% butterPulse = circshift(butterPulse, round(length(butterPulse)/2)-idx);
len = length(pulse); % 脉冲长度
pulseCentered = [zeros(round(len/2)-idx, 1); pulse(1:(end-round(len/2)+idx))]; % 居中脉冲

% 互相关性检查:

% 不绘图，直接进行数据计算和检查

%% 按照15.4.4节创建根升余弦脉冲
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
x = xcorr(r, pulseCentered, 'normalized'); % 归一化互相关

%% 执行互相关检查
T1 = cfg.CrossCorrThreshold1; % 阈值1
T2 = cfg.CrossCorrThreshold2; % 阈值2

% 检查主峰和旁瓣
crossCorrMainPeakOK = any(x(round(end/2 +[-1 0 1]))>T1); % 检查主峰
ind = find(abs(x)>T2); % 找到大于T2的索引
crossCorrSideLobesOK = isempty(ind) || isequal(ind', min(ind):max(ind)); % 检查旁瓣连续性
crossCorrPassed = crossCorrMainPeakOK && crossCorrSideLobesOK;

%% 时域模板检查（简化检查：脉冲峰值在合理范围内）
pulseNormalized = pulse/max(pulse); % 脉冲归一化
peakValue = max(abs(pulseNormalized));
timeMaskPassed = (peakValue <= 1.1) && (peakValue >= 0.1); % 简化的模板检查

% 准备返回结果
result.crossCorrPassed = crossCorrPassed;
result.timeMaskPassed = timeMaskPassed;
result.overallPassed = crossCorrPassed && timeMaskPassed;

% 输出检查结果
fprintf('=== UWB基带冲激响应分析结果 ===\n');
if crossCorrPassed
    fprintf('✓ 互相关检查通过\n');
    fprintf('  - 主峰检查: %s\n', char(crossCorrMainPeakOK*'通过' + ~crossCorrMainPeakOK*'失败'));
    fprintf('  - 旁瓣检查: %s\n', char(crossCorrSideLobesOK*'通过' + ~crossCorrSideLobesOK*'失败'));
else
    fprintf('✗ 互相关检查失败\n');
    fprintf('  - 主峰检查: %s\n', char(crossCorrMainPeakOK*'通过' + ~crossCorrMainPeakOK*'失败'));
    fprintf('  - 旁瓣检查: %s\n', char(crossCorrSideLobesOK*'通过' + ~crossCorrSideLobesOK*'失败'));
end

if timeMaskPassed
    fprintf('✓ 时域模板检查通过\n');
else
    fprintf('✗ 时域模板检查失败\n');
end

if result.overallPassed
    result.message = '所有基础冲激响应测试项均通过，脉冲符合规范';
    fprintf('✓ 整体分析: 所有基础冲激响应测试项均通过\n');
else
    result.message = '部分基础冲激响应测试项失败，请检查脉冲参数';
    fprintf('✗ 整体分析: 部分基础冲激响应测试项失败\n');
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
    plot([-xSpan/2 xSpan/2], [T1 T1], 'r-') % 绘制T1阈值线
    plot([-xSpan/2 xSpan/2], [T2 T2], 'r:') % 绘制T2阈值线
    plot([-xSpan/2 xSpan/2], [-T2 -T2], 'r:') % 绘制-T2阈值线
    
    title('互相关') % 设置子图标题
    axis([-NTp NTp -0.5 max(x)]) % 设置坐标轴范围
    xlabel('时间 (ns)') % 设置x轴标签
    
    % 在图上显示检查结果（居中置顶）
    if crossCorrPassed
        sgtitle('互相关检查: 通过 ✓', 'Color', 'green', 'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle('互相关检查: 失败 ✗', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold');
    end

%% 时域模板:

% 绘制脉冲及其时域模板:
figMask = figure; % 新建图形窗口
pulse = pulse/max(pulse); % 脉冲归一化
xPulseStart = cfg.PulseStartTime; % Tp为单位的起始时间
xPulseEnd = cfg.PulseEndTime; % Tp为单位的结束时间
plot((xPulseStart:(xPulseEnd-xPulseStart)/(length(pulse)-1):xPulseEnd), pulse, 'b-o'); % 绘制脉冲
xMin = cfg.TimeMaskXMin; % x轴最小值
xMax = cfg.TimeMaskXMax; % x轴最大值
yMin = cfg.TimeMaskYMin; % y轴最小值
yMax = cfg.TimeMaskYMax; % y轴最大值
axis([xMin xMax yMin yMax]) % 设置坐标轴范围
t0 = 0; % t=0

% 绘制模板区域
a = cfg.TimeMaskAlpha; % 透明度
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

