function [btype,analog,errStr,msgobj] = uwb_iirchk(Wn,varargin)
%IIRCHK  Parameter checking for BUTTER, CHEBY1, CHEBY2, and ELLIP.
%   [btype,analog,errStr] = iirchk(Wn,varargin) returns the
%   filter type btype (1=lowpass, 2=bandpss, 3=highpass, 4=bandstop)
%   and analog flag analog (0=digital, 1=analog) given the edge
%   frequency Wn (either a one or two element vector) and the
%   optional arguments in varargin.  The variable arguments are
%   either empty, a one element cell, or a two element cell.
%
%   errStr is empty if no errors are detected; otherwise it contains
%   the error message.  If errStr is not empty, btype and analog
%   are invalid.

%   Copyright 1988-2019 The MathWorks, Inc.
%#codegen
errStr = '';
msgobj = [];

if coder.target('MATLAB')
    % Define defaults:
    analog = 0; % 0=digital, 1=analog
    btype = 1;  % 1=lowpass, 2=bandpss, 3=highpass, 4=bandstop

    if length(Wn)==1
        btype = 1;
    elseif length(Wn)==2
        btype = 2;
    else
        msgobj = message('signal:iirchk:MustBeOneOrTwoElementVector','Wn');
        errStr = getString(msgobj);
        return
    end

    if length(varargin)>2
        msgobj = message('signal:iirchk:TooManyInputArguments');
        errStr = getString(msgobj);
        return
    end

    % Interpret and strip off trailing 's' or 'z' argument:
    if ~isempty(varargin)
        switch lower(varargin{end})
          case 's'
            analog = 1;
            varargin(end) = [];
          case 'z'
            analog = 0;
            varargin(end) = [];
          otherwise
            if length(varargin) > 1
                msgobj = message('signal:iirchk:BadAnalogFlag','z','s');
                errStr = getString(msgobj);
                return
            end
        end
    end

    % Check for correct Wn limits
    if ~analog
        if any(Wn<=0) | any(Wn>=1)
            msgobj = message('signal:iirchk:FreqsMustBeWithinUnitInterval');
            errStr = getString(msgobj);
            return
        end
    else
        if any(Wn<=0)
            msgobj = message('signal:iirchk:FreqsMustBePositive');
            errStr = getString(msgobj);
            return
        end
    end

    % At this point, varargin will either be empty, or contain a single
    % band type flag.

    if length(varargin)==1   % Interpret filter type argument:
        switch lower(varargin{1})
          case 'low'
            btype = 1;
          case 'bandpass'
            btype = 2;
          case 'high'
            btype = 3;
          case 'stop'
            btype = 4;
          otherwise
            if nargin == 2
                msgobj = message('signal:iirchk:BadOptionString', ...
                                 'high','stop','low','bandpass','z','s');
                errStr = getString(msgobj);
            else  % nargin == 3
                msgobj = message('signal:iirchk:BadFilterType', ...
                                 'high','stop','low','bandpass');
                errStr = getString(msgobj);
            end
            return
        end
        switch btype
          case 1
            if length(Wn)~=1
                msgobj = message('signal:iirchk:BadOptionLength','low','Wn',1);
                errStr = getString(msgobj);
                return
            end
          case 2
            if length(Wn)~=2
                msgobj = message('signal:iirchk:BadOptionLength','bandpass','Wn',2);
                errStr = getString(msgobj);
                return
            end
          case 3
            if length(Wn)~=1
                msgobj = message('signal:iirchk:BadOptionLength','high','Wn',1);
                errStr = getString(msgobj);
                return
            end
          case 4
            if length(Wn)~=2
                msgobj = message('signal:iirchk:BadOptionLength','stop','Wn',2);
                errStr = getString(msgobj);
                return
            end
        end
    end

else
    coder.internal.errorIf(length(Wn) ~= 1 && length(Wn) ~= 2,'signal:iirchk:MustBeOneOrTwoElementVector','Wn')
    coder.internal.errorIf(length(varargin) > 2,'signal:iirchk:TooManyInputArguments');
    for i = 1:length(varargin)
        validateattributes(varargin{i},{'char'},{'vector'},'iirchk','',i+1);
    end
    % Interpret trailing 's' or 'z' argument:
    if ~isempty(varargin)
        analogFlag = lower(varargin{end});
        isAnalogFlag = (strcmp(analogFlag,'s') || strcmp(analogFlag,'z'));
        coder.internal.errorIf(~isAnalogFlag && length(varargin) == 2,'signal:iirchk:BadAnalogFlag','z','s');
        if strcmp(analogFlag,'s')
            analog = true;
            idx = 1;
        elseif strcmp(analogFlag,'z')
            analog = false;
            idx = 1;
        else
            analog = false;
            idx = 0;
        end
    else
        idx = 0;
        analog = false;
    end
    isFilterSpecified = logical(length(varargin) - idx);
    % analogFlag = 's' or 'z', filtType = 'high','low','bandpass','stop'
    % 1.iirchk(wn)                     ->  length(varargin) = 0, idx = 0 => isFilterSpecified = 0;
    % 2.iirchk(wn,analogFlag)          ->  length(varargin) = 1, idx = 1 => isFilterSpecified = 0
    % 3.iirchk(wn,filtType)            ->  length(varargin) = 1, idx = 0 => isFilterSpecified = 1;
    % 4.iirchk(wn,filtType,analogFlag) ->  length(varargin) = 2, idx = 1 => isFilterSpecified = 1;
    % 5.iirchk(wn,incorrectInput)      ->  length(varargin) = 1, idx = 0 => isFilterSpecified = 1;
    % For case 5 we throw error from line 167

    % Check for correct Wn limits
    coder.internal.errorIf(~analog && (any(Wn <= 0,'all') || any(Wn >= 1,'all')),'signal:iirchk:FreqsMustBeWithinUnitInterval');
    coder.internal.errorIf(analog && any(Wn <= 0,'all'),'signal:iirchk:FreqsMustBePositive');

    if isFilterSpecified    % Interpret filter type argument:
        filter = lower(varargin{1});
        badFilter = ~(strcmp(filter,'low') || strcmp(filter,'high') || strcmp(filter,'bandpass') || strcmp(filter,'stop'));
        coder.internal.errorIf(~isAnalogFlag && badFilter && nargin == 2,'signal:iirchk:BadOptionString','high','stop','low','bandpass','z','s');
        coder.internal.errorIf(badFilter && nargin == 3,'signal:iirchk:BadFilterType','high','stop','low','bandpass');
        if ~badFilter
            if strcmp(filter,'low')
                btype = 1;
            elseif strcmp(filter,'bandpass')
                btype = 2;
            elseif strcmp(filter,'high')
                btype = 3;
            else % stop
                btype = 4;
            end
        else
            btype = 0; % This is only for coder inference
        end
        coder.internal.errorIf(btype == 1 && length(Wn) ~= 1,'signal:iirchk:BadOptionLength','low','Wn',1);
        coder.internal.errorIf(btype == 2 && length(Wn) ~= 2,'signal:iirchk:BadOptionLength','bandpass','Wn',2);
        coder.internal.errorIf(btype == 3 && length(Wn) ~= 1,'signal:iirchk:BadOptionLength','high','Wn',1);
        coder.internal.errorIf(btype == 4 && length(Wn) ~= 2,'signal:iirchk:BadOptionLength','stop','Wn',2);
    else
        if length(Wn) == 1
            btype = 1;
        else
            btype = 2;
        end
    end
end

% LocalWords:  CHEBY btype Wn bandpss wn
