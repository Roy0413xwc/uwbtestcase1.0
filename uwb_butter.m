function varargout = uwb_butter(n, Wn, varargin)
%BUTTER Butterworth digital and analog filter design.
%   [B,A] = BUTTER(N,Wn) designs an Nth order lowpass digital
%   Butterworth filter and returns the filter coefficients in length
%   N+1 vectors B (numerator) and A (denominator). The coefficients
%   are listed in descending powers of z. The cutoff frequency
%   Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to
%   half the sample rate.
%
%   If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an
%   order 2N bandpass filter with passband  W1 < W < W2.
%   [B,A] = BUTTER(N,Wn,'high') designs a highpass filter.
%   [B,A] = BUTTER(N,Wn,'low') designs a lowpass filter.
%   [B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
%
%   When used with three left-hand arguments, as in
%   [Z,P,K] = BUTTER(...), the zeros and poles are returned in
%   length N column vectors Z and P, and the gain in scalar K.
%
%   When used with four left-hand arguments, as in
%   [A,B,C,D] = BUTTER(...), state-space matrices are returned.
%
%   BUTTER(N,Wn,'s'), BUTTER(N,Wn,'high','s') and BUTTER(N,Wn,'stop','s')
%   design analog Butterworth filters.  In this case, Wn is in [rad/s]
%   and it can be greater than 1.0.
%
%   % Example 1:
%   %   For data sampled at 1000 Hz, design a 9th-order highpass
%   %   Butterworth filter with cutoff frequency of 300Hz.
%
%   Wn = 300/500;                   % Normalized cutoff frequency
%   [z,p,k] = butter(9,Wn,'high');  % Butterworth filter
%   [b1,a1] = zp2ctf(z,p,k);        % Convert to CTF form
%   filterAnalyzer(b1,a1);          % Plot magnitude response
%
%   % Example 2:
%   %   Design a 4th-order butterworth band-pass filter which passes
%   %   frequencies between 0.15 and 0.3.
%
%   [b2,a2]=butter(2,[.15,.3]);      % Bandpass digital filter design
%   h = filterAnalyzer(b2,a2);       % Visualize filter
%
%   See also BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ,
%   FILTER, DESIGNFILT.

%   Copyright 1998-2023 The MathWorks, Inc.

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

%#codegen

narginchk(2,4);
if coder.target('MATLAB')
   [varargout{1:nargout}] = butterImpl(n,Wn,varargin{:});
else
   allConst = coder.internal.isConst(n) && coder.internal.isConst(Wn);
   for ii = 1:length(varargin)
      allConst = allConst && coder.internal.isConst(varargin{ii});
   end
   if allConst && coder.internal.isCompiled
      [varargout{1:nargout}] = coder.const(@feval,'butter',n,Wn,varargin{:});
   else
      [varargout{1:nargout}] = butterImpl(n,Wn,varargin{:});
   end
end
end

function varargout = butterImpl(n,Wn,varargin)
inputArgs = cell(1,length(varargin));
if nargin > 2
   [inputArgs{:}] = convertStringsToChars(varargin{:});
else
   inputArgs = varargin;
end
validateattributes(n,{'numeric'},{'scalar','real','integer','positive'},'butter','N');
validateattributes(Wn,{'numeric'},{'vector','real','finite','nonempty'},'butter','Wn');

[btype,analog,~,msgobj] = uwb_iirchk(Wn,inputArgs{:});
if ~isempty(msgobj)
   coder.internal.error(msgobj.Identifier,msgobj.Arguments{:});
end
% Cast to enforce precision rules
n1 = double(n(1));
coder.internal.errorIf(n1 > 500,'signal:butter:InvalidRange')
% Cast to enforce precision rules
Wn = double(Wn);
% step 1: get analog, pre-warped frequencies
fs = 2;
if ~analog
   u = 2*fs*tan(pi*Wn/fs);
else
   u = Wn;
end

% step 2: Get N-th order Butterworth analog lowpass prototype
[zs,ps,ks] = buttap(n1);
% Transform to state-space
[a,b,c,d] = zp2ss(zs,ps,ks);
% step 3: Transform to the desired filter
if length(Wn) == 1
   % step 3a: convert to low-pass prototype estimate
   Wn1 = u(1);
   Bw = []; %#ok<NASGU>
   % step 3b: Transform to lowpass or high pass filter of desired cutoff
   % frequency
   if btype == 1           % Lowpass
      [ad,bd,cd,dd] = lp2lp(a,b,c,d,Wn1);
   else % btype == 3       % Highpass
      [ad,bd,cd,dd] = lp2hp(a,b,c,d,Wn1);
   end
else % length(Wn) is 2
   % step 3a: convert to low-pass prototype estimate
   Bw = u(2) - u(1);      % center frequency
   Wn1 = sqrt(u(1)*u(2));
   % step 3b: Transform to bandpass or bandstop filter of desired center
   % frequency and bandwidth
   if btype == 2           % Bandpass
      [ad,bd,cd,dd] = lp2bp(a,b,c,d,Wn1,Bw);
   else % btype == 4       % Bandstop
      [ad,bd,cd,dd] = lp2bs(a,b,c,d,Wn1,Bw);
   end
end
% step 4: Use Bilinear transformation to find discrete equivalent:
if ~analog
   [ad,bd,cd,dd] = bilinear(ad,bd,cd,dd,fs);
end

if nargout == 4 % Outputs are in state space form
   varargout{1} = ad;          % A
   varargout{2} = bd;          % B
   varargout{3} = cd;          % C
   varargout{4} = dd;          % D
else
   p = eig(ad);
   [z,k] = buttzeros(btype,n1,Wn1,analog,p+0i);
   if nargout == 3         % Transform to zero-pole-gain form
      varargout{1} = z;
      varargout{2} = p;
      varargout{3} = k;
   else
      den = real(poly(p));
      num = [zeros(1,length(p)-length(z),'like',den)  k*real(poly(z))];
      varargout{1} = num;
      varargout{2} = den;
   end
end
end


function [z,k] = buttzeros(btype,n,Wn,analog,p)
% This internal function computes the zeros and gain of the ZPK
% representation. Wn is scalar (sqrt(Wn(1)*Wn(2)) for bandpass/stop).
if analog
   % for lowpass and bandpass, don't include zeros at +Inf or -Inf
   switch btype
      case 1  % lowpass: H(0)=1
         z = zeros(0,1,'like',p);
         k = Wn^n;  % prod(-p) = Wn^n
      case 2  % bandpass: H(1i*Wn) = 1
         z = zeros(n,1,'like',p);
         k = real(prod(1i*Wn-p)/(1i*Wn)^n);
      case 3  % highpass: H(Inf) = 1
         z = zeros(n,1,'like',p);
         k = 1;
      case 4  % bandstop: H(0) = 1
         z = 1i*Wn*((-1).^(0:2*n-1)');
         k = 1;  % prod(p) = prod(z) = Wn^(2n)
      otherwise
         coder.internal.error('signal:iirchk:BadFilterType','high','stop','low','bandpass');
   end
else
   Wn = 2*atan2(Wn,4);
   switch btype
      case 1  % lowpass: H(1)=1
         z = -ones(n,1,'like',p);
         k = real(prod(1-p))/2^n;
      case 2  % bandpass: H(z) = 1 for z=exp(1i*sqrt(Wn(1)*Wn(2)))
         z = [ones(n,1,'like',p); -ones(n,1,'like',p)];
         zWn = exp(1i*Wn);
         k = real(prod(zWn-p)/prod(zWn-z));
      case 3  % highpass: H(-1) = 1
         z = ones(n,1,'like',p);
         k = real(prod(1+p))/2^n;
      case 4  % bandstop: H(1) = 1
         z = exp(1i*Wn*( (-1).^(0:2*n-1)' ));
         k = real(prod(1-p)/prod(1-z));
      otherwise
         coder.internal.error('signal:iirchk:BadFilterType','high','stop','low','bandpass');
   end
end
% Note: codegen complains when z set to both real and complex values above
if ~any(imag(z))
   z = real(z);
end
end


% LocalWords:  Butterworth Wn th butterworth CHEBY DESIGNFILT btype infzero
% LocalWords:  iirchk Shure Krauss Burrus
