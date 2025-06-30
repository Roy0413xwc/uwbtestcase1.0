% lrwpanHRPConfig 属性:
% Channel                 - 信道编号
% Mode                    - 工作模式
% MeanPRF                 - 平均脉冲重复频率(单位:MHz)
% DataRate                - 有效载荷数据速率(单位:Mbps)
% PHRDataRate             - PHR数据速率(单位:Mbps)
% SamplesPerPulse         - 每个巴特沃斯脉冲的采样点数
% STSPacketConfiguration  - 控制STS在数据包中的位置
% NumSTSSegments          - STS段的数量
% STSSegmentLength        - 有效STS段长度(以512码片为倍数)
% ExtraSTSGapLength       - 额外STS间隔长度(以4码片为倍数)
% ExtraSTSGapIndex        - 额外STS间隔值的索引
% CodeIndex               - 从表15-6、15-7和15-7a中选择的SYNC码索引
% PreambleMeanPRF         - 前导码的平均PRF(脉冲重复频率，单位:MHz)
% PreambleDuration        - 扩频前导SYNC码的重复次数
% SFDNumber               - 帧起始定界符(SFD)选择索引
% Ranging                 - 标记PHY帧是否用于测距的标志
% ConstraintLength        - 卷积编码的约束长度标志
% PSDULength              - PHY服务数据单元长度(单位:字节)
% SampleRate              - 波形采样率
%
%
% PeakPRF                 - 峰值脉冲重复频率
% BurstsPerSymbol         - 每个符号的突发数(见表15-3)
% NumHopBursts            - 每个符号的候选突发数(见表15-3)
% ChipsPerBurst           - 每个突发的码片数(见表15-3)
% ChipsPerSymbol          - 每个符号的码片数(见表15-3)
% ConvolutionalCoding     - 有效载荷卷积编码器激活标志
% PreambleCodeLength      - 前导码的符号长度
% PreambleSpreadingFactor - 增量函数扩频码长度
%
% 测试项参数:
% ButterworthOrder        - Butterworth滤波器阶数
% ButterworthCutoff       - Butterworth滤波器截止频率(单位:Hz)
% PulseWidth              - 脉冲宽度(单位:ns)
% PulseWidthMultiplier    - 脉冲宽度倍数
% RRCBeta                 - RRC滚降系数
% TimeMaskXMin            - 时域模板x轴最小值(Tp单位)
% TimeMaskXMax            - 时域模板x轴最大值(Tp单位)
% TimeMaskYMin            - 时域模板y轴最小值
% TimeMaskYMax            - 时域模板y轴最大值
% TimeMaskAlpha           - 模板透明度
% PulseStartTime          - 脉冲起始时间(Tp单位)
% PulseEndTime            - 脉冲结束时间(Tp单位)
% CrossCorrThreshold1     - 互相关主峰检查阈值
% CrossCorrThreshold2     - 互相关旁瓣检查阈值
% PSDPeakPRF              - PSD模板峰值脉冲重复频率(单位:Hz)
% PSDFc065Factor          - PSD模板0.65倍峰值频率因子
% PSDDrop065              - PSD模板0.65倍频率处的下降(单位:dB)
% PSDFc08Factor           - PSD模板0.8倍峰值频率因子
% PSDDrop08               - PSD模板0.8倍频率处的下降(单位:dB)
% PSDRB                   - PSD分析分辨带宽(单位:Hz)
% PSDYLimits              - PSD分析Y轴范围(单位:dB)

function cfg = uwb_paramsconfig()

 cfg= struct(...
    ...%%————————————————————发送端参数——————————————————————————
    'Channel', 0,  ...                   % 通道号，范围[0:3 5:6 8:10 12:14]
    'Mode', 'HPRF',   ...                % 工作模式：'HPRF'、'BPRF'或'802.15.4a'
    'MeanPRF', 249.6,  ...               % 平均脉冲重复频率（MHz）
    'DataRate', 0.85,    ...             % 有效载荷数据速率（Mbps）
    'PHRDataRate', 0.85,  ...            % PHY头数据速率（Mbps）
    'SamplesPerPulse', 4, ...            % 每个脉冲的采样点数
    'STSPacketConfiguration', 1, ...     % STS在数据包中的位置（0-3）
    'NumSTSSegments', 1,   ...           % STS段数量（HPRF模式）
    'STSSegmentLength', 64,   ...        % STS段码片长度（512码片为单位）
    'ExtraSTSGapLength', 0,  ...         % 额外STS间隔长度（4码片倍数）
    'ExtraSTSGapIndex', 0,  ...          % 额外STS间隔索引
    'CodeIndex', 25,    ...              % SYNC码索引
    'PreambleMeanPRF', 16.1,  ...        % 前导码平均PRF（MHz）
    'PreambleDuration', 64,  ...         % 前导码重复次数
    'SFDNumber', 0,  ...                 % 帧起始定界符索引
    'Ranging', false,   ...              % 是否用于测距
    'ConstraintLength', 3,  ...          % 卷积编码约束长度
    'PSDULength', 127,  ...              % PHY服务数据单元长度（字节）
    ...% 依赖参数（计算值）
    'PeakPRF', [],   ...                 % 峰值脉冲重复频率（MHz）
    'BurstsPerSymbol', [],  ...          % 每个符号的突发数
    'NumHopBursts', [],  ...             % 每个符号的候选突发数
    'ChipsPerBurst', [],  ...            % 每个突发的码片数
    'ChipsPerSymbol', [],  ...           % 每个符号的码片数
    'ConvolutionalCoding', [], ...       % 卷积编码启用标志
    'PreambleCodeLength', [], ...        % 前导码符号长度
    'PreambleSpreadingFactor', [], ...   % 扩频因子
    'SampleRate', [],  ...               % 输出采样率（Hz）
    'MeanPRFNum', [],  ...               % 内部使用的MeanPRF映射值
    'DataRateNum', [],   ...             % 内部使用的DataRate映射值
    ...%%————————————————————测试项参数——————————————————————————
    ...% Butterworth脉冲参数
    'ButterworthOrder', 4,  ...          % Butterworth滤波器阶数
    'ButterworthCutoff', 500e6,  ...     % Butterworth滤波器截止频率 (Hz)
    'PulseWidth', 2,  ...                % 脉冲宽度 (ns)
    'PulseWidthMultiplier', 8,  ...      % 脉冲宽度倍数
    ...% RRC脉冲参数
    'RRCBeta', 0.5,  ...                 % RRC滚降系数
    ...% 时域模板参数
    'TimeMaskXMin', -3,  ...             % 时域模板x轴最小值 (Tp单位)
    'TimeMaskXMax', 9,  ...              % 时域模板x轴最大值 (Tp单位)
    'TimeMaskYMin', -0.8,  ...           % 时域模板y轴最小值
    'TimeMaskYMax', 1.1,  ...            % 时域模板y轴最大值
    'TimeMaskAlpha', 0.75,  ...          % 模板透明度
    'PulseStartTime', -1.25,  ...        % 脉冲起始时间 (Tp单位)
    'PulseEndTime', 8,  ...              % 脉冲结束时间 (Tp单位)
    ...% 互相关检查参数
    'CrossCorrThreshold1', 0.8,  ...     % 主峰检查阈值
    'CrossCorrThreshold2', 0.3,  ...     % 旁瓣检查阈值
    ...% 传输PSD遮罩参数
    'PSDPeakPRF', 499.2e6,  ...          % 峰值脉冲重复频率 (Hz)
    'PSDFc065Factor', 0.65,  ...         % 0.65倍峰值频率因子
    'PSDDrop065', -10,  ...              % 0.65倍频率处的下降 (dB)
    'PSDFc08Factor', 0.8,  ...           % 0.8倍峰值频率因子
    'PSDDrop08', -18,  ...               % 0.8倍频率处的下降 (dB)
    'PSDRB', 1e6,  ...                   % 分辨带宽 (Hz)
    'PSDYLimits', [-350, -40]  ...       % Y轴范围 (dB)
    ...
    ...%%————————————————————接收端参数——————————————————————————
);
end


