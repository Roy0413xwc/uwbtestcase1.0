function cfg = uwb_init_Params(cfg)
% 初始化依赖参数，根据输入参数推导其他只读字段

%% ————————————————————————发送端初始化————————————————————

  % 1. PeakPRF
  peak = 499.2;
  if strcmp(cfg.Mode, 'HPRF') && cfg.MeanPRF == 124.8
    peak = peak / 2;
  end
  cfg.PeakPRF = peak;

  % 2. MeanPRFNum
  if strcmp(cfg.Mode, 'BPRF')
    cfg.MeanPRFNum = 62.4;
  else
    cfg.MeanPRFNum = cfg.MeanPRF;
  end

  % 3. DataRateNum
  if strcmp(cfg.Mode, 'BPRF')
    cfg.DataRateNum = 6810;
  else
    cfg.DataRateNum = cfg.DataRate * 1e3; % 单位 kbps
  end

  % 4. BurstsPerSymbol
  switch cfg.MeanPRFNum
    case 15.6
      cfg.BurstsPerSymbol = 32;
    case 3.9
      cfg.BurstsPerSymbol = 128;
    otherwise % 默认62.4
      cfg.BurstsPerSymbol = 8;
  end

  % 5. NumHopBursts
  cfg.NumHopBursts = cfg.BurstsPerSymbol / 4;

  % 6. ChipsPerBurst
  prf = cfg.MeanPRFNum;
  rate = cfg.DataRateNum;
  if prf == 15.6
    if rate == 110
      cfg.ChipsPerBurst = 128;
    elseif rate == 850
      cfg.ChipsPerBurst = 16;
    elseif rate == 6810
      cfg.ChipsPerBurst = [16, 2];
    else
      cfg.ChipsPerBurst = [16, 1];
    end
  elseif prf == 3.9
    if rate == 110
      cfg.ChipsPerBurst = 32;
    elseif rate == 850
      cfg.ChipsPerBurst = 4;
    elseif rate == 1700
      cfg.ChipsPerBurst = [4, 2];
    else
      cfg.ChipsPerBurst = [4, 1];
    end
  else % 62.4 MHz
    if rate == 110
      cfg.ChipsPerBurst = 512;
    elseif rate == 850
      cfg.ChipsPerBurst = 64;
    elseif rate == 6810
      if cfg.PHRDataRate == 0.85
        cfg.ChipsPerBurst = [64, 8];
      else
        cfg.ChipsPerBurst = 8;
      end
    else
      cfg.ChipsPerBurst = [64, 2];
    end
  end

  % 7. ChipsPerSymbol
  if strcmp(cfg.Mode, 'HPRF')
    if cfg.ConstraintLength == 3
      cfg.ChipsPerSymbol = [16 8] * 249.6 / cfg.MeanPRFNum;
    else
      cfg.ChipsPerSymbol = 8 * 249.6 / cfg.MeanPRFNum;
    end
  else
    cfg.ChipsPerSymbol = cfg.ChipsPerBurst * cfg.BurstsPerSymbol;
  end
  % 8. ConvolutionalCoding
  cfg.ConvolutionalCoding = ~((cfg.MeanPRFNum == 3.9 && rate == 6810) || ...
                                    (cfg.MeanPRFNum == 15.6 && rate == 27240));

  % 9. PreambleCodeLength
  if cfg.CodeIndex <= 8
    cfg.PreambleCodeLength = 31;
  elseif cfg.CodeIndex <= 24
    cfg.PreambleCodeLength = 127;
  else
    cfg.PreambleCodeLength = 91;
  end

  % 10. PreambleSpreadingFactor
  if cfg.PreambleCodeLength > 31
    cfg.PreambleSpreadingFactor = 4;
  else
    if cfg.PreambleMeanPRF == 16.1
      cfg.PreambleSpreadingFactor = 16;
    else
      cfg.PreambleSpreadingFactor = 64;
    end
  end

  % 11. SampleRate (Hz)
  cfg.SampleRate = cfg.SamplesPerPulse * cfg.PeakPRF * 1e6;

  %% ————————————————————————接收端初始化————————————————————



end
