function SFD = uwb_getSFD(cfg)
% lrwpan.internal.getSFD 获取适用的帧起始定界符(SFD)序列
% SFD = lrpwan.internal.getSFD(CFG) 根据输入的<a
% href="matlab:help('lrwpanHRPConfig')">lrwpanHRPConfig</a>配置对象CFG中指定的
% 操作模式(HPRF、BPRF或802.15.4a)和SFD索引或数据速率，返回对应的SFD序列。
% 输出序列SFD的长度为4、8、16、32或64个符号。

%   Copyright 2021-2023 The MathWorks, Inc.

%#codegen

if ~strcmp(cfg.Mode, '802.15.4a')
  % 见翻译标准表16-11
  switch cfg.SFDNumber
    case 0
      SFD = [0 +1 0 -1 +1 0 0 -1]; % legacy
    case 1
      SFD = [-1 -1 +1 -1];
    case 2
      SFD = [-1 -1 -1 +1 -1 -1 +1 -1];
    case 3
      SFD = [-1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1 -1 +1 -1];
    otherwise % case 4
      SFD = [-1 -1 -1 -1 -1 -1 -1 +1 -1 -1 +1 -1 -1 +1 -1 +1 -1 +1 ...
-1 -1 -1 +1 +1 -1 -1 -1 +1 -1 +1 +1 -1 -1];
  end  
  
else % 15.4a
  % 见标准16.2.6.3
  if cfg.DataRateNum ~= 110
    SFD = [0 +1 0 -1 +1 0 0 -1]; % short SFD
  else % long SFD
    SFD = [0 +1 0 -1 +1 0 0 -1 0 +1 0 -1 +1 0 0 -1 -1 0 0 +1 0 -1 0 +1 0 +1 0 ...
      0 0 -1 0 -1 0 -1 0 0 +1 0 -1 -1 0 -1 +1 0 0 0 0 +1 +1 0 0 -1 -1 -1 +1 -1 ...
      +1 +1 0 0 0 0 +1 +1];
  end
end
