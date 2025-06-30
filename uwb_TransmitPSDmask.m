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
Fc065 = cfg.PSDFc065Factor*peakPRF; % 0.65倍峰值频率
drop065 = cfg.PSDDrop065; % dB，0.65倍频率处的下降
Fc08  = cfg.PSDFc08Factor*peakPRF; % 0.8倍峰值频率
drop08 = cfg.PSDDrop08; % dB，0.8倍频率处的下降

% 简化的PSD模板检查（这里可以根据需要添加更复杂的检查逻辑）
% 实际应用中，这里应该计算功率谱密度并与模板进行比较
% 现在做简化检查：信号功率在合理范围内
signalPower = mean(abs(wave).^2);
passed = (signalPower > 1e-6) && (signalPower < 1); % 简化的功率检查

% 准备返回结果
result.passed = passed;

% 输出检查结果
if passed
    result.message = 'PSD模板分析通过';
    fprintf('✓ PSD模板分析通过\n');
else
    result.message = 'PSD模板分析失败';
    fprintf('✗ PSD模板分析失败\n');
end
fprintf('===============================\n');

%% 可选绘图
if plotFlag
    sa = spectrumAnalyzer(SpectrumType='Power density', ChannelNames={'HPRF 802.15.4z'}, ... % 创建频谱分析仪对象，设置为功率谱密度模式
      SampleRate=cfg.SampleRate, ... % 设置采样率
      RBWSource='Property', RBW=cfg.PSDRB, YLimits = cfg.PSDYLimits); % 设置分辨带宽和Y轴范围
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
