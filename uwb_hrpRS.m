function out = uwb_hrpRS(in, doEncode)
%  lrwpan.internal.hrpRS HRP物理层的Reed-Solomon编码/解码
%  OUT = lrpwan.internal.hrpRS(IN, DOENCODE) 根据IEEE Std 802.15.4™‐2020标准第15.3.3.2节执行Reed-Solomon编码或解码。
%  IN是任意长度的二进制列向量。DOENCODE是一个标志，当其值分别为true或false时，决定执行编码或解码操作。
%  当DOENCODE为true时，OUT是编码输出；当DOENCODE为false时，OUT是解码输出。

%   Copyright 2021-2023 The MathWorks, Inc.

% As per Sec. 15.3.3.2 in IEEE Std 802.15.4™‐2020

%#codegen

persistent rsEnc rsDec

N = 63;%编码后的码块长度  55->63
K = 55;%信息符号数
M = 6; %每个符号的比特率
genPoly = 'x8 + 55x7 + 61x6 + 37x5 + 48x4 + 47x3 + 20x2 + 6x1 + 22';
primPoly = '1 + x + x6';

out = []; % init output

if isempty(rsEnc) && doEncode
  rsEnc = comm.RSEncoder(N, K, genPoly, 'BitInput', false,  ... % do not use BitInput because comm.RS*coder use left-msb
    'PrimitivePolynomialSource', "Property", 'PrimitivePolynomial', primPoly);
end

if isempty(rsDec) && ~doEncode
  rsDec = comm.RSDecoder(N, K, genPoly, 'BitInput', false,  ...
    'PrimitivePolynomialSource', "Property", 'PrimitivePolynomial', primPoly);
end

blockSize = 330;                      % 330 for encoding  将其扩展为330个零（哑元）比特
if ~doEncode
  blockSize = blockSize + M*(63-55);  % 378 for decoding
end

% Process PSDU in blocks of 330 bits
for blockIdx = 1:ceil(length(in)/blockSize)
  thisBlock = in(1+blockSize*(blockIdx-1) : min(end, blockSize*blockIdx));
  I = length(thisBlock);

  % a) Addition of dummy bits
  inPadded = [zeros(blockSize-I, 1); thisBlock];
  
  % b) Bit-to-symbol conversion, with right-msb
  msbFirst = false;
  inPaddedInt = bit2int(inPadded, M, msbFirst);
  
  % c) Encoding/Decoding
  if doEncode
    tmp = rsEnc(inPaddedInt(1:(blockSize/M)));  % index to help codegen
  else
    tmp = rsDec(inPaddedInt(1:(blockSize/M)));
  end
  
  % d) Symbol to bit conversion
  outPaddedBits = int2bit(tmp, M, msbFirst); % 
  
  % e) Removal of dummy bits, concatenation with previous outputs
  out = [out; outPaddedBits((1 + max(0, blockSize-I)) :end)]; %#ok<AGROW>
end

