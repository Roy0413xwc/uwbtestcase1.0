function uwb_validate_Params(cfg)
% 验证结构体 params 中的参数是否满足 IEEE 802.15.4z 要求
% 若检测到无效配置，函数将抛出错误

%% ————————————————————————发送端参数校验————————————————————

  % Mode vs MeanPRF
  if (strcmp(cfg.Mode, 'HPRF') && any(cfg.MeanPRFNum == [3.9, 15.6, 62.4])) || ...
     (strcmp(cfg.Mode, '802.15.4a') && any(cfg.MeanPRFNum == [124.8, 249.6]))
    error('Invalid combination of Mode and MeanPRFNum.');
  end

  % MeanPRF vs DataRate
  if (cfg.MeanPRFNum == 3.9 && cfg.DataRateNum == 27240) || ...
     (cfg.MeanPRFNum > 3.9 && cfg.DataRateNum == 1700)
    error('Invalid combination of MeanPRFNum and DataRateNum.');
  end

  % CodeIndex vs MeanPRF
  if (cfg.CodeIndex <= 8 && cfg.MeanPRFNum >= 62.4) || ...
     (cfg.CodeIndex > 8 && cfg.CodeIndex < 25 && cfg.MeanPRFNum ~= 62.4) || ...
     (cfg.CodeIndex >= 25 && strcmp(cfg.Mode, '802.15.4a'))
    error('Invalid combination of CodeIndex and MeanPRFNum.');
  end

  % CodeIndex vs Channel
  switch cfg.CodeIndex
    case {1, 2}
      allowedChannels = [0, 1, 8, 12];
    case {3, 4}
      allowedChannels = [2, 5, 9, 13];
    case {5, 6}
      allowedChannels = [3, 6, 10, 14];
    otherwise
      allowedChannels = setdiff(0:15, [4, 7, 11, 15]);
  end

  if ~ismember(cfg.Channel, allowedChannels)
    error('Channel %d not allowed for CodeIndex %d.', cfg.Channel, cfg.CodeIndex);
  end

  % Legacy mode check
  if cfg.MeanPRFNum == 62.4 && cfg.DataRateNum ~= 6810
    if cfg.CodeIndex > 25
      error('Invalid CodeIndex for legacy 62.4 MHz configuration.');
    end
  end

  % PreambleDuration vs Mode & CodeIndex
  if (strcmp(cfg.Mode, 'HPRF') && ismember(cfg.PreambleDuration, [1024, 4096])) || ...
     (~strcmp(cfg.Mode, 'HPRF') && ismember(cfg.PreambleDuration, [24, 32, 48, 96, 128, 256])) || ...
     (~strcmp(cfg.Mode, 'HPRF') && cfg.PreambleDuration == 4096 && ...
      cfg.CodeIndex <= 8 && cfg.PreambleMeanPRF == 4.03)
    error('Invalid PreambleDuration for current configuration.');
  end

  % Mode vs SFDNumber
  if strcmp(cfg.Mode, 'BPRF') && ~ismember(cfg.SFDNumber, [0, 2])
    error('Invalid SFDNumber for BPRF mode.');
  end

  % PSDULength
  if (~strcmp(cfg.Mode, 'HPRF') && cfg.PSDULength > 127) || ...
     (strcmp(cfg.Mode, 'HPRF') && cfg.STSPacketConfiguration == 2 && ...
      cfg.ExtraSTSGapLength > 0 && cfg.PSDULength > 1023)
    error('Invalid PSDULength for current configuration.');
  end


   %% ————————————————————————接收端参数校验————————————————————



end
