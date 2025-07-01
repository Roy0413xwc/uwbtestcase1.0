function result = uwb_TransmitPSDmask(wave, cfg, plotFlag)
% uwb_TransmitPSDmask - UWB传输功率谱密度模板分析
% 这是一个内部示例辅助函数，随时可能更改。
%
% 输入:
%   wave - 波形数据
%   cfg - 配置参数结构体
%   plotFlag - 可选参数，是否绘制图形 (默认为false)
%
% 输出:
%   result - 结构体，包含检查结果:
%            .passed: PSD模板检查是否通过 (logical)
%            .message: 检查结果描述

%   版权所有 2021-2023 The MathWorks, Inc

if nargin < 3
    plotFlag = false; % 默认不绘图
end

% 计算PSD模板参数
peakPRF = cfg.PSDPeakPRF; % 峰值脉冲重复频率

% PSD模板固定参数（符合IEEE 802.15.4a标准）
PSD_FC065_FACTOR = 0.65;  % 0.65倍峰值频率因子
PSD_DROP065 = -10;        % 0.65倍频率处的下降 (dB)
PSD_FC08_FACTOR = 0.8;    % 0.8倍峰值频率因子  
PSD_DROP08 = -18;         % 0.8倍频率处的下降 (dB)
PSD_RBW = 1e6;            % 分辨带宽 (Hz)

Fc065 = PSD_FC065_FACTOR*peakPRF; % 0.65倍峰值频率
drop065 = PSD_DROP065; % dB，0.65倍频率处的下降
Fc08  = PSD_FC08_FACTOR*peakPRF; % 0.8倍峰值频率
drop08 = PSD_DROP08; % dB，0.8倍频率处的下降

% 使用Welch谱估计法计算功率谱密度
% Welch方法：对长度为N的数据分段，允许每段有50%重叠，用窗函数平滑处理
N = length(wave);
nfft = 2^nextpow2(N/8);  % FFT点数，通常取N/8的2次幂
window = hann(nfft);     % 汉宁窗函数
noverlap = floor(nfft/2); % 50%重叠
fs = cfg.SampleRate;     % 采样频率

% 计算Welch功率谱密度
[psd, freq] = pwelch(wave, window, noverlap, nfft, fs, 'psd');

% 转换为dBm/MHz单位 (假设阻抗为50欧姆)
psd_dBm_MHz = 10*log10(psd * 1000 * 1e6); % 转换为dBm/MHz

% 找到峰值频率和功率
[psd_peak_dBm, peak_idx] = max(psd_dBm_MHz);
peak_freq = freq(peak_idx);

% 计算模板频率点
Fc065_actual = Fc065; % 0.65倍峰值频率
Fc08_actual = Fc08;   % 0.8倍峰值频率

% 在PSD中找到对应频率点的功率值
[~, idx_065] = min(abs(freq - Fc065_actual));
[~, idx_08] = min(abs(freq - Fc08_actual));

% 获取模板频率点处的实际功率值（相对于峰值）
power_at_065 = psd_dBm_MHz(idx_065) - psd_peak_dBm; % 相对功率
power_at_08 = psd_dBm_MHz(idx_08) - psd_peak_dBm;   % 相对功率

% PSD模板合规性检查
tolerance = 2; % 容差2dB
passed_065 = power_at_065 <= (drop065 + tolerance); % 0.65倍频率处检查
passed_08 = power_at_08 <= (drop08 + tolerance);    % 0.8倍频率处检查

% 整体检查结果
passed = passed_065 && passed_08;

% 准备返回结果
result.passed = passed;
result.psd = psd;           % 功率谱密度
result.freq = freq;         % 频率向量
result.psd_peak_dBm = psd_peak_dBm;  % 峰值功率
result.peak_freq = peak_freq;        % 峰值频率
result.power_at_065 = power_at_065;  % 0.65倍频率处相对功率
result.power_at_08 = power_at_08;    % 0.8倍频率处相对功率

% 输出检查结果
if passed
    result.message = sprintf(['PSD模板分析通过\n' ...
                             '峰值功率: %.2f dBm/MHz @ %.2f MHz\n' ...
                             '0.65倍频率处: %.2f dB (要求: <= %.2f dB)\n' ...
                             '0.8倍频率处: %.2f dB (要求: <= %.2f dB)'], ...
                             psd_peak_dBm, peak_freq/1e6, ...
                             power_at_065, drop065, ...
                             power_at_08, drop08);
    fprintf('✓ PSD模板分析通过\n');
    fprintf('  峰值功率: %.2f dBm/MHz @ %.2f MHz\n', psd_peak_dBm, peak_freq/1e6);
    fprintf('  0.65倍频率处: %.2f dB (要求: <= %.2f dB) %s\n', ...
            power_at_065, drop065, char(10004)); % ✓
    fprintf('  0.8倍频率处: %.2f dB (要求: <= %.2f dB) %s\n', ...
            power_at_08, drop08, char(10004)); % ✓
else
    result.message = sprintf(['PSD模板分析失败\n' ...
                             '峰值功率: %.2f dBm/MHz @ %.2f MHz\n' ...
                             '0.65倍频率处: %.2f dB (要求: <= %.2f dB) %s\n' ...
                             '0.8倍频率处: %.2f dB (要求: <= %.2f dB) %s'], ...
                             psd_peak_dBm, peak_freq/1e6, ...
                             power_at_065, drop065, char(10060*passed_065 + 10004*(~passed_065)), ...
                             power_at_08, drop08, char(10060*passed_08 + 10004*(~passed_08)));
    fprintf('✗ PSD模板分析失败\n');
    fprintf('  峰值功率: %.2f dBm/MHz @ %.2f MHz\n', psd_peak_dBm, peak_freq/1e6);
    fprintf('  0.65倍频率处: %.2f dB (要求: <= %.2f dB) %s\n', ...
            power_at_065, drop065, char(10060*passed_065 + 10004*(~passed_065)));
    fprintf('  0.8倍频率处: %.2f dB (要求: <= %.2f dB) %s\n', ...
            power_at_08, drop08, char(10060*passed_08 + 10004*(~passed_08)));
end
fprintf('===============================\n');

%% Welch功率谱密度图形展示
if plotFlag
    % 创建Welch PSD图形
    figure('Name', 'Welch功率谱密度分析', 'Position', [100, 100, 800, 600]);
    
    % 绘制Welch计算的PSD
    plot(freq/1e6, psd_dBm_MHz, 'b-', 'LineWidth', 1.5);
    hold on;
    
    % 标记峰值点
    plot(peak_freq/1e6, psd_peak_dBm, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    % 标记关键频率点
    plot(Fc065_actual/1e6, psd_dBm_MHz(idx_065), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(Fc08_actual/1e6, psd_dBm_MHz(idx_08), 'ms', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
    
    % 绘制PSD模板线
    % 构建模板频率和功率向量
    template_freq = [-inf, -Fc08_actual, -Fc065_actual, 0, Fc065_actual, Fc08_actual, inf] / 1e6; % MHz
    template_power = [drop08, drop08, drop065, 0, drop065, drop08, drop08] + psd_peak_dBm; % dBm/MHz
    
    % 只绘制正频率部分的模板
    pos_idx = template_freq >= 0;
    plot(template_freq(pos_idx), template_power(pos_idx), 'r--', 'LineWidth', 2, 'DisplayName', 'PSD模板');
    
    % 绘制负频率部分的模板（镜像）
    neg_idx = template_freq <= 0;
    plot(-template_freq(neg_idx), template_power(neg_idx), 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    % 设置图形属性
    xlabel('频率 (MHz)');
    ylabel('功率谱密度 (dBm/MHz)');
    title('UWB信号Welch功率谱密度分析');
    grid on;
    xlim([0, max(freq)/1e6]);
    
    % 添加图例
    legend({'Welch PSD', sprintf('峰值 (%.2f MHz)', peak_freq/1e6), ...
            sprintf('0.65f_c (%.2f MHz)', Fc065_actual/1e6), ...
            sprintf('0.8f_c (%.2f MHz)', Fc08_actual/1e6), ...
            'PSD模板'}, 'Location', 'best');
    
    % 添加文本标注
    text(0.02, 0.95, sprintf('峰值功率: %.2f dBm/MHz', psd_peak_dBm), ...
         'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
    text(0.02, 0.90, sprintf('0.65倍频率处: %.2f dB', power_at_065), ...
         'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
    text(0.02, 0.85, sprintf('0.8倍频率处: %.2f dB', power_at_08), ...
         'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
    text(0.02, 0.80, sprintf('检查结果: %s', char(passed*10004 + ~passed*10060)), ...
         'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    
    hold off;
    
    % 创建第二个图形：显示Welch方法的详细信息
    figure('Name', 'Welch方法参数信息', 'Position', [920, 100, 600, 400]);
    
    % 创建子图1：时域信号
    subplot(2,1,1);
    t = (0:length(wave)-1) / fs * 1e6; % 时间轴（微秒）
    plot(t, real(wave), 'b-', 'LineWidth', 1);
    xlabel('时间 (μs)');
    ylabel('幅度');
    title('原始UWB信号');
    grid on;
    
    % 创建子图2：Welch方法参数展示
    subplot(2,1,2);
    % 绘制窗函数
    window_demo = window(1:min(length(window), 500)); % 限制显示长度
    plot(window_demo, 'r-', 'LineWidth', 1.5);
    xlabel('样本点');
    ylabel('窗函数值');
    title(sprintf('Welch参数: 窗长=%d, 重叠=%d%%, FFT点数=%d', ...
                  length(window), round(noverlap/length(window)*100), nfft));
    grid on;
    
    % 添加参数信息文本
    info_text = {
        sprintf('信号长度: %d 样本', N),
        sprintf('采样频率: %.2f MHz', fs/1e6),
        sprintf('窗函数: Hann窗'),
        sprintf('窗长度: %d 样本', length(window)),
        sprintf('重叠长度: %d 样本 (%.1f%%)', noverlap, noverlap/length(window)*100),
        sprintf('FFT点数: %d', nfft),
        sprintf('频率分辨率: %.2f kHz', fs/nfft/1000)
    };
    
    % 在图形上添加信息文本
    annotation('textbox', [0.02, 0.02, 0.4, 0.4], 'String', info_text, ...
               'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
end

%% 可选绘图（原有的频谱分析仪）
if plotFlag
    sa = spectrumAnalyzer(SpectrumType='Power density', ChannelNames={'HPRF 802.15.4z'}, ... % 创建频谱分析仪对象，设置为功率谱密度模式
      SampleRate=cfg.SampleRate, ... % 设置采样率
      RBWSource='Property', RBW=PSD_RBW, YLimits = cfg.PSDYLimits); % 设置分辨带宽和Y轴范围
    sa(wave); % 显示输入信号的频谱

    spectralMask = SpectralMaskSpecification; % 创建频谱模板对象
    spectralMask.EnabledMasks='Upper'; % 只启用上模板
    spectralMask.ReferenceLevel = 'Spectrum peak'; % 参考电平为频谱峰值
    spectralMask.UpperMask = [-inf drop08; -Fc08 drop08; -Fc08 drop065; -Fc065 drop065; -Fc065 0; % 设置上模板的各个点
                              Fc065 0; Fc065 drop065; Fc08 drop065; Fc08 drop08; inf drop08];
    sa.SpectralMask = spectralMask; % 应用频谱模板到分析仪

    assignin('caller', 'spectrumAnalyzer', sa); % 将频谱分析仪对象分配到工作区
end

end
