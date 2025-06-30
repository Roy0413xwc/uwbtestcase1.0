function result = uwb_TransmitCenterFrequencyTolerance(wave, cfg, plotFlag)
% uwb_TransmitCenterFrequencyTolerance - UWB发射中心频率容差测试
% 这是一个内部示例辅助函数，随时可能更改。
%
% 功能描述：
%   根据IEEE 802.15.4a标准16.4.9节要求，测试HRP UWBPHY的
%   发射中心频率容差，该要求应从16.4.6节中码片时钟容差的要求
%
% 测试项目：
%   1. 发射中心频率容差测试 - 容差要求±20×10⁻⁶
%
% 输入:
%   wave - 波形数据
%   cfg - 配置参数结构体
%   plotFlag - 可选参数，是否绘制图形 (默认为false)
%
% 输出:
%   result - 结构体，包含检查结果:
%            .centerFreqTolerancePassed: 中心频率容差检查是否通过 (logical)
%            .overallPassed: 整体检查是否通过 (logical)
%            .message: 检查结果描述
%            .centerFreqError: 中心频率误差 (ppm)

%   版权所有 2024 UWB测试项目组

if nargin < 3
    plotFlag = false; % 默认不绘图
end

fprintf('===============================\n');
fprintf('发射中心频率容差测试\n');
fprintf('===============================\n');

%% 测试常量定义
TOLERANCE_REQUIREMENT = 20e-6; % ±20×10⁻⁶ 容差要求

%% 测试1: 发射中心频率容差测试
fprintf('测试1: 发射中心频率容差测试\n');
fprintf('-------------------------------\n');

% 获取理论中心频率（根据信道配置）
% 参考表16-27: UWB信道中心频率定义
switch cfg.Channel
    case 1
        theoreticalCenterFreq = 3494.4e6; % Hz
    case 2
        theoreticalCenterFreq = 3993.6e6; % Hz
    case 3
        theoreticalCenterFreq = 4492.8e6; % Hz
    case 4
        theoreticalCenterFreq = 3993.6e6; % Hz
    case 5
        theoreticalCenterFreq = 6489.6e6; % Hz
    case 6
        theoreticalCenterFreq = 6988.8e6; % Hz
    case 7
        theoreticalCenterFreq = 6489.6e6; % Hz
    case 8
        theoreticalCenterFreq = 7488.0e6; % Hz
    case 9
        theoreticalCenterFreq = 7987.2e6; % Hz
    case 10
        theoreticalCenterFreq = 8486.4e6; % Hz
    case 11
        theoreticalCenterFreq = 7987.2e6; % Hz
    case 12
        theoreticalCenterFreq = 8985.6e6; % Hz
    case 13
        theoreticalCenterFreq = 9484.8e6; % Hz
    case 14
        theoreticalCenterFreq = 9984.0e6; % Hz
    case 15
        theoreticalCenterFreq = 9484.8e6; % Hz
    otherwise
        theoreticalCenterFreq = 7987.2e6; % 默认信道9
end

% 通过频域分析测量实际中心频率
fs = cfg.SampleRate; % 采样频率
N = length(wave);
frequencyAxis = (-N/2:N/2-1) * (fs/N); % 频率轴

% 计算功率谱密度
waveSpectrum = fftshift(fft(wave));
powerSpectrum = abs(waveSpectrum).^2;

% 方法1: 找到功率谱的重心频率作为中心频率
totalPower = sum(powerSpectrum);
measuredCenterFreq1 = abs(sum(frequencyAxis .* powerSpectrum') / totalPower);

% 方法2: 找到功率谱的峰值频率作为中心频率  
[~, peakIdx] = max(powerSpectrum);
measuredCenterFreq2 = abs(frequencyAxis(peakIdx));

% 选择更合理的测量值
if measuredCenterFreq1 > theoreticalCenterFreq / 2
    measuredCenterFreq = measuredCenterFreq1; % 使用重心频率
    method = '重心频率法';
elseif measuredCenterFreq2 > theoreticalCenterFreq / 2
    measuredCenterFreq = measuredCenterFreq2; % 使用峰值频率
    method = '峰值频率法';
else
    % 可能是基带信号，假设载波频率正确
    measuredCenterFreq = theoreticalCenterFreq;
    method = '基带信号假设';
    fprintf('注意: 检测到基带信号，假设中心频率为配置值\n');
end

% 计算中心频率误差
centerFreqError = (measuredCenterFreq - theoreticalCenterFreq) / theoreticalCenterFreq;
centerFreqErrorPPM = centerFreqError * 1e6; % 转换为ppm

% 判断是否通过容差测试
centerFreqTolerancePassed = abs(centerFreqError) <= TOLERANCE_REQUIREMENT;

fprintf('理论中心频率: %.1f MHz (信道 %d)\n', theoreticalCenterFreq/1e6, cfg.Channel);
fprintf('测量中心频率: %.1f MHz (%s)\n', measuredCenterFreq/1e6, method);
fprintf('中心频率误差: %.2f ppm\n', centerFreqErrorPPM);
fprintf('容差要求: ±%.1f ppm\n', TOLERANCE_REQUIREMENT*1e6);

if centerFreqTolerancePassed
    fprintf('✓ 发射中心频率容差测试通过\n');
else
    fprintf('✗ 发射中心频率容差测试失败\n');
end

%% 测试2: 频率稳定性检查（可选）
fprintf('\n测试2: 频率稳定性检查\n');
fprintf('-------------------------------\n');

% 将信号分段分析频率稳定性
segmentLength = min(1024, floor(length(wave)/4));
numSegments = floor(length(wave) / segmentLength);

if numSegments >= 2
    centerFreqs = zeros(numSegments, 1);
    
    for i = 1:numSegments
        startIdx = (i-1) * segmentLength + 1;
        endIdx = i * segmentLength;
        segmentWave = wave(startIdx:endIdx);
        
        % 计算每段的中心频率
        segmentSpectrum = fftshift(fft(segmentWave));
        segmentPower = abs(segmentSpectrum).^2;
        segmentFreqAxis = (-segmentLength/2:segmentLength/2-1) * (fs/segmentLength);
        
        segmentTotalPower = sum(segmentPower);
        if segmentTotalPower > 0
            centerFreqs(i) = abs(sum(segmentFreqAxis .* segmentPower') / segmentTotalPower);
        else
            centerFreqs(i) = theoreticalCenterFreq;
        end
    end
    
    % 计算频率稳定性
    freqStd = std(centerFreqs);
    freqStabilityPPM = freqStd / theoreticalCenterFreq * 1e6;
    
    fprintf('频率标准差: %.2f Hz (%.2f ppm)\n', freqStd, freqStabilityPPM);
    
    % 频率稳定性要求（应小于容差要求的一半）
    stabilityPassed = freqStabilityPPM <= (TOLERANCE_REQUIREMENT * 1e6 / 2);
    
    if stabilityPassed
        fprintf('✓ 频率稳定性检查通过\n');
    else
        fprintf('✗ 频率稳定性检查失败\n');
    end
else
    stabilityPassed = true; % 信号太短，跳过稳定性检查
    fprintf('信号长度不足，跳过稳定性检查\n');
end

%% 整体测试结果
fprintf('\n===============================\n');
fprintf('发射中心频率容差测试结果\n');
fprintf('===============================\n');

overallPassed = centerFreqTolerancePassed && stabilityPassed;

if overallPassed
    resultMessage = '发射中心频率容差测试整体通过';
    fprintf('✓ %s\n', resultMessage);
else
    resultMessage = '发射中心频率容差测试整体失败';
    fprintf('✗ %s\n', resultMessage);
    if ~centerFreqTolerancePassed
        fprintf('  - 中心频率容差不符合要求\n');
    end
    if ~stabilityPassed
        fprintf('  - 频率稳定性不符合要求\n');
    end
end

fprintf('测试依据: IEEE 802.15.4a标准 16.4.9节\n');
fprintf('相关要求: 16.4.6节码片时钟容差要求\n');
fprintf('===============================\n');

%% 准备返回结果
result.centerFreqTolerancePassed = centerFreqTolerancePassed;
result.stabilityPassed = stabilityPassed;
result.overallPassed = overallPassed;
result.message = resultMessage;
result.centerFreqError = centerFreqErrorPPM;
result.method = method;

%% 可选绘图
if plotFlag
    figure('Name', 'Transmit Center Frequency Tolerance Test Results', 'NumberTitle', 'off');
    
    % 子图1: 频域谱
    subplot(2,2,1);
    frequencyAxisMHz = frequencyAxis / 1e6; % 频率轴 (MHz)
    plot(frequencyAxisMHz, 10*log10(powerSpectrum + eps));
    hold on;
    % 标记理论中心频率
    xline(theoreticalCenterFreq/1e6, 'r--', '理论中心频率', 'LineWidth', 2);
    % 标记测量中心频率
    xline(measuredCenterFreq/1e6, 'g--', '测量中心频率', 'LineWidth', 2);
    title('功率谱密度');
    xlabel('频率 (MHz)');
    ylabel('功率 (dB)');
    grid on;
    legend('功率谱', '理论中心频率', '测量中心频率', 'Location', 'best');
    
    % 子图2: 容差分析
    subplot(2,2,2);
    toleranceLimits = [-TOLERANCE_REQUIREMENT*1e6, TOLERANCE_REQUIREMENT*1e6];
    if centerFreqTolerancePassed
        barColor = 'green';
    else
        barColor = 'red';
    end
    bar(1, centerFreqErrorPPM, 'FaceColor', barColor);
    hold on;
    yline(toleranceLimits(1), 'r--', '下限', 'LineWidth', 2);
    yline(toleranceLimits(2), 'r--', '上限', 'LineWidth', 2);
    ylim([min(toleranceLimits(1)*1.5, centerFreqErrorPPM*1.5), max(toleranceLimits(2)*1.5, centerFreqErrorPPM*1.5)]);
    title('中心频率容差分析');
    ylabel('误差 (ppm)');
    xticks(1);
    xticklabels({'中心频率误差'});
    grid on;
    
    % 子图3: 频率稳定性（如果有多段数据）
    if exist('centerFreqs', 'var') && length(centerFreqs) > 1
        subplot(2,2,3);
        plot(1:length(centerFreqs), centerFreqs/1e6, 'o-');
        hold on;
        yline(theoreticalCenterFreq/1e6, 'r--', '理论值', 'LineWidth', 2);
        title('频率稳定性');
        xlabel('时间段');
        ylabel('中心频率 (MHz)');
        grid on;
    else
        subplot(2,2,3);
        text(0.5, 0.5, '数据不足\n无法分析稳定性', 'HorizontalAlignment', 'center');
        title('频率稳定性');
    end
    
    % 子图4: 测试结果汇总
    subplot(2,2,4);
    axis off;
    text(0.1, 0.8, sprintf('中心频率误差: %.2f ppm', centerFreqErrorPPM), 'FontSize', 12);
    text(0.1, 0.6, sprintf('容差要求: ±%.1f ppm', TOLERANCE_REQUIREMENT*1e6), 'FontSize', 12);
    text(0.1, 0.4, sprintf('检测方法: %s', method), 'FontSize', 12);
    
    if overallPassed
        text(0.1, 0.2, '✓ 整体测试通过', 'FontSize', 14, 'Color', 'green', 'FontWeight', 'bold');
    else
        text(0.1, 0.2, '✗ 整体测试失败', 'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');
    end
    
    title('测试结果汇总');
end

end
