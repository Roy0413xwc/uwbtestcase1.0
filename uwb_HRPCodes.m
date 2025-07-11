function code = uwb_HRPCodes(codeIndex)
%  lrwpan.internal.HRPCodes 获取特定前导码(SYNC)序列
%  CODE = lrpwan.internal.HRPCodes(CODEINDEX) 根据标准16.2.6.2节的规定，返回由CODEINDEX索引的特定前导码(SYNC)序列。
%  当CODEINDEX在1到8之间时，CODE为31符号序列，见表格16-7。
%  当CODEINDEX在9到24之间时，CODE为127符号序列，见表格16-8。
%  当CODEINDEX在25到32之间时，CODE为91符号序列，见表格16-9。

%   Copyright 2021-2023 The MathWorks, Inc.

%#codegen

persistent Codes

if isempty(Codes)

  Codes = { %% Length 31 - Table 15-6
    [-1 0 0 0 0 +1 0 -1 0 +1 +1 +1 0 +1 -1 0 0 0 +1 -1 +1 +1 +1 0 0 -1 +1 0 -1 0 0]; % 1
    [0 +1 0 +1 -1 0 +1 0 +1 0 0 0 -1 +1 +1 0 -1 +1 -1 -1 -1 0 0 +1 0 0 +1 +1 0 0 0]; % 2
    [-1 +1 0 +1 +1 0 0 0 -1 +1 -1 +1 +1 0 0 +1 +1 0 +1 0 0 -1 0 0 0 0 -1 0 +1 0 -1]; % 3
    [0 0 0 0 +1 -1 0 0 -1 0 0 -1 +1 +1 +1 +1 0 +1 -1 +1 0 0 0 +1 0 -1 0 +1 +1 0 -1]; % 4
    [-1 0 +1 -1 0 0 +1 +1 +1 -1 +1 0 0 0 -1 +1 0 +1 +1 +1 0 -1 0 +1 0 0 0 0 -1 0 0]; % 5
    [+1 +1 0 0 +1 0 0 -1 -1 -1 +1 -1 0 +1 +1 -1 0 0 0 +1 0 +1 0 -1 +1 0 +1 0 0 0 0]; % 6
    [+1 0 0 0 0 +1 -1 0 +1 0 +1 0 0 +1 0 0 0 +1 0 +1 +1 -1 -1 -1 0 -1 +1 0 0 -1 +1]; % 7
    [0 +1 0 0 -1 0 -1 0 +1 +1 0 0 0 0 -1 -1 +1 0 0 -1 +1 0 +1 +1 -1 +1 +1 0 +1 0 0]; % 8
    
    
    %% Length 127 - Table 15-7
    [+1 0 0 +1 0 0 0 -1 0 -1 -1 0 0 -1 -1 +1 0 +1 0 +1 0 0 -1 +1 -1 +1 +1 0 +1 0 0 ...
    0 0 +1 +1 -1 0 0 0 +1 0 0 -1 0 0 -1 -1 0 -1 +1 0 +1 0 -1 -1 0 -1 +1 +1 +1 0 +1 ...
    +1 0 0 0 +1 -1 0 +1 0 0 -1 0 +1 +1 -1 0 +1 +1 +1 0 0 -1 +1 0 0 +1 0 +1 0 -1 0 ...
    +1 +1 -1 +1 -1 -1 +1 0 0 0 0 0 0 +1 0 0 0 0 0 -1 +1 0 0 0 0 -1 0 -1 0 0 0 -1 -1 +1]; %9
    
    [+1 +1 0 0 +1 0 -1 +1 0 0 +1 0 0 +1 0 0 0 0 0 0 -1 0 0 0 -1 0 0 -1 -1 0 0 0 -1 0 +1 ...
    -1 +1 0 -1 0 +1 -1 0 -1 +1 0 0 0 0 0 +1 -1 0 0 +1 +1 0 -1 0 +1 0 0 -1 -1 +1 0 0 +1 ...
    +1 -1 +1 0 +1 -1 0 +1 0 0 0 0 -1 0 -1 0 -1 0 -1 +1 +1 -1 +1 0 +1 0 0 +1 0 +1 0 0 0 ...
    -1 +1 0 +1 +1 +1 0 0 0 -1 -1 -1 -1 +1 +1 +1 0 0 0 0 +1 +1 +1 0 -1 -1];             %10
    
    [-1 +1 -1 0 0 0 0 +1 0 0 -1 -1 0 0 0 0 0 -1 0 +1 0 +1 0 +1 -1 0 +1 0 0 +1 0 0 +1 ...
    0 -1 0 0 -1 +1 +1 +1 0 0 +1 0 0 0 -1 +1 0 +1 0 -1 0 0 0 0 +1 +1 +1 +1 +1 -1 +1 0 ...
    +1 -1 -1 0 +1 -1 0 +1 +1 -1 -1 0 -1 0 0 0 +1 0 -1 +1 0 0 +1 0 +1 -1 -1 -1 -1 0 0 ...
    0 -1 0 0 0 0 0 0 -1 +1 0 0 +1 -1 0 +1 +1 0 0 0 +1 +1 -1 0 0 +1 +1 -1 0 -1 0];      % 11
    
    [-1 +1 0 +1 +1 0 0 0 0 0 0 -1 0 +1 0 -1 +1 0 -1 -1 -1 +1 -1 +1 +1 0 0 -1 +1 0 ...
    +1 +1 0 +1 0 +1 0 +1 0 0 0 -1 0 0 -1 0 0 -1 +1 0 0 +1 -1 +1 +1 0 0 0 -1 +1 -1 0 ...
    -1 +1 +1 0 -1 0 +1 +1 +1 +1 0 -1 0 0 -1 0 +1 +1 0 0 +1 0 +1 0 0 +1 +1 -1 0 0 ...
    +1 0 0 0 +1 -1 0 0 0 -1 0 -1 -1 +1 0 0 0 0 -1 0 0 0 0 -1 -1 0 +1 0 0 0 0 0 +1 -1 -1]; %12
    
    [+1 0 0 0 -1 -1 0 0 0 0 -1 -1 +1 +1 0 -1 +1 +1 +1 +1 0 -1 0 +1 +1 0 +1 0 -1 0 0 ...
    -1 +1 0 +1 +1 0 0 +1 +1 -1 0 +1 +1 0 +1 -1 +1 0 -1 0 0 +1 0 0 -1 0 -1 -1 0 0 0  ...
    -1 +1 -1 0 0 +1 0 0 0 0 -1 0 +1 +1 -1 0 0 0 0 0 +1 -1 0 -1 0 0 0 0 0 0 -1 0 0 -1 ...
    +1 -1 +1 +1 -1 +1 0 0 0 -1 0 +1 0 +1 0 +1 +1 +1 -1 0 0 -1 -1 0 0 +1 0 +1 0 0 0];     %13
    
    [+1 0 0 0 +1 +1 0 -1 0 +1 0 -1 0 0 +1 -1 0 -1 +1 0 -1 0 0 +1 0 +1 0 0 0 0 +1 0 ...
    +1 -1 0 0 0 0 +1 +1 0 0 +1 0 +1 +1 +1 +1 +1 -1 +1 0 -1 0 +1 -1 0 -1 -1 +1 0 +1 ...
    +1 -1 -1 0 0 0 -1 -1 -1 0 +1 0 0 0 +1 0 +1 0 -1 +1 -1 0 0 0 0 0 0 +1 -1 +1 -1 ...
    0 -1 -1 0 0 +1 +1 0 0 0 -1 0 0 +1 0 0 +1 +1 -1 0 0 -1 -1 +1 +1 -1 0 0 -1 0 0 0 0 0]; %14
    
    [0 +1 -1 0 0 +1 0 -1 0 0 0 -1 +1 +1 0 0 0 0 -1 -1 -1 +1 +1 0 0 0 +1 0 +1 -1 0 -1 ...
    +1 0 0 -1 +1 0 0 0 -1 -1 0 -1 0 0 -1 -1 0 -1 -1 +1 +1 +1 -1 +1 0 -1 +1 +1 0 0 +1 ...
    -1 +1 +1 0 +1 0 0 0 0 0 +1 0 -1 0 +1 +1 +1 -1 0 0 +1 0 0 +1 0 0 0 -1 0 0 0 0 +1 ...
    0 0 -1 -1 +1 0 +1 +1 0 +1 0 +1 0 -1 0 0 -1 0 -1 +1 -1 0 +1 0 +1 +1 0 0 0 0 0];     % 15
    
    [+1 +1 0 0 0 0 +1 0 0 0 +1 0 0 +1 -1 -1 0 +1 -1 +1 +1 0 -1 0 0 0 -1 -1 0 0 +1 -1 ...
    0 +1 0 0 +1 +1 0 0 0 +1 +1 +1 0 0 +1 0 +1 0 -1 0 -1 +1 -1 0 -1 0 +1 0 0 +1 0 0 +1 ...
    0 +1 +1 -1 -1 -1 -1 +1 0 0 +1 +1 -1 -1 +1 0 +1 -1 0 -1 -1 +1 0 0 0 0 0 0 -1 0 -1 ...
    0 0 0 0 -1 +1 0 -1 -1 0 0 +1 0 0 0 0 0 +1 -1 +1 +1 0 0 0 -1 0 -1 +1 0 +1 0];       % 16
    
    [+1 -1 -1 0 0 0 -1 0 -1 0 0 0 0 +1 -1 0 0 0 0 0 +1 0 0 0 0 0 0 +1 -1 -1 +1 -1 +1 ...
    +1 0 -1 0 +1 0 +1 0 0 +1 -1 0 0 +1 +1 +1 0 -1 +1 +1 0 -1 0 0 +1 0 -1 +1 0 0 0 +1 ...
    +1 0 +1 +1 +1 -1 0 -1 -1 0 +1 0 +1 -1 0 -1 -1 0 0 -1 0 0 +1 0 0 0 -1 +1 +1 0 0 0 ...
    0 +1 0 +1 +1 -1 +1 -1 0 0 +1 0 +1 0 +1 -1 -1 0 0 -1 -1 0 -1 0 0 0 +1 0 0 +1];      % 17
    
    [-1 -1 0 +1 +1 +1 0 0 0 0 +1 +1 +1 -1 -1 -1 -1 0 0 0 +1 +1 +1 0 +1 -1 0 0 0 +1 ...
    0 +1 0 0 +1 0 +1 -1 +1 +1 -1 0 -1 0 -1 0 -1 0 0 0 0 +1 0 -1 +1 0 +1 -1 +1 +1 0 ...
    0 +1 -1 -1 0 0 +1 0 -1 0 +1 +1 0 0 -1 +1 0 0 0 0 0 +1 -1 0 -1 +1 0 -1 0 +1 -1 +1 ...
    0 -1 0 0 0 -1 -1 0 0 -1 0 0 0 -1 0 0 0 0 0 0 +1 0 0 +1 0 0 +1 -1 0 +1 0 0 +1 +1];   % 18
    
    [-1 0 -1 +1 +1 0 0 -1 +1 +1 0 0 0 +1 +1 0 -1 +1 0 0 +1 -1 0 0 0 0 0 0 -1 0 0 0 ...
    -1 -1 -1 -1 +1 0 +1 0 0 +1 -1 0 +1 0 0 0 -1 0 -1 -1 +1 +1 0 -1 +1 0 -1 -1 +1 0 ...
    +1 -1 +1 +1 +1 +1 +1 0 0 0 0 -1 0 +1 0 +1 -1 0 0 0 +1 0 0 +1 +1 +1 -1 0 0 -1 ...
    0 +1 0 0 +1 0 0 +1 0 -1 +1 0 +1 0 +1 0 -1 0 0 0 0 0 -1 -1 0 0 +1 0 0 0 0 -1 +1 -1 0]; % 19
    
    [-1 -1 +1 0 0 0 0 0 +1 0 -1 -1 0 0 0 0 -1 0 0 0 0 +1 -1 -1 0 -1 0 0 0 -1 +1 0 0 0 ...
    +1 0 0 -1 +1 +1 0 0 +1 0 +1 0 0 +1 +1 0 -1 0 0 -1 0 +1 +1 +1 +1 0 -1 0 +1 +1 -1 0 ...
    -1 +1 -1 0 0 0 +1 +1 -1 +1 0 0 +1 -1 0 0 -1 0 0 -1 0 0 0 +1 0 +1 0 +1 0 +1 +1 0 +1 ...
    -1 0 0 +1 +1 -1 +1 -1 -1 -1 0 +1 -1 0 +1 0 -1 0 0 0 0 0 0 +1 +1 0 +1 -1];             % 20
    
    [+1 0 +1 0 0 -1 -1 0 0 -1 +1 +1 +1 0 +1 0 +1 0 -1 0 0 0 +1 -1 +1 +1 -1 +1 -1 0 0 ...
    -1 0 0 0 0 0 0 -1 0 -1 +1 0 0 0 0 0 -1 +1 +1 0 -1 0 0 0 0 +1 0 0 -1 +1 -1 0 0 0 ...
    -1 -1 0 -1 0 0 +1 0 0 -1 0 +1 -1 +1 0 +1 +1 0 -1 +1 +1 0 0 +1 +1 0 +1 -1 0 0 -1 ...
    0 +1 0 +1 +1 0 -1 0 +1 +1 +1 +1 -1 0 +1 +1 -1 -1 0 0 0 0 -1 -1 0 0 0 +1 0 0 0];       % 21
    
    [0 -1 0 0 -1 +1 +1 -1 -1 0 0 -1 +1 +1 0 0 +1 0 0 -1 0 0 0 +1 +1 0 0 -1 -1 0 -1 +1 ...
    -1 +1 0 0 0 0 0 0 -1 +1 -1 0 +1 0 +1 0 0 0 +1 0 -1 -1 -1 0 0 0 -1 -1 +1 +1 0 +1 -1 ...
    -1 0 -1 +1 0 -1 0 +1 -1 +1 +1 +1 +1 +1 0 +1 0 0 +1 +1 0 0 0 0 -1 +1 0 +1 0 0 0 0 +1 ...
    0 +1 0 0 -1 0 +1 -1 0 -1 +1 0 0 -1 0 +1 0 -1 0 +1 +1 0 0 0 +1 0 0 0 0];               % 22
    
    [0 0 0 +1 +1 0 +1 0 -1 +1 -1 0 -1 0 0 -1 0 +1 0 +1 0 +1 +1 0 +1 -1 -1 0 0 +1 0 0 0 ...
    0 -1 0 0 0 +1 0 0 +1 0 0 -1 +1 +1 +1 0 -1 0 +1 0 0 0 0 0 +1 0 +1 +1 -1 +1 0 0 +1 +1 ...
    -1 0 +1 -1 +1 +1 +1 -1 -1 0 -1 -1 0 0 -1 0 -1 -1 0 0 0 +1 -1 0 0 +1 -1 0 -1 +1 0 +1 ...
    0 0 0 +1 +1 -1 -1 -1 0 0 0 0 +1 +1 -1 0 0 0 -1 0 +1 0 0 -1 +1 0 0 0];                 % 23
    
    [+1 0 +1 -1 0 -1 0 0 0 +1 +1 -1 +1 0 0 0 0 0 +1 0 0 -1 -1 0 +1 -1 0 0 0 0 -1 0 -1 0 ...
    0 0 0 0 0 +1 -1 -1 0 -1 +1 0 +1 -1 -1 +1 +1 0 0 +1 -1 -1 -1 -1 +1 +1 0 +1 0 0 +1 0 0 ...
    +1 0 -1 0 -1 +1 -1 0 -1 0 +1 0 +1 0 0 +1 +1 +1 0 0 0 +1 +1 0 0 +1 0 -1 +1 0 0 -1 -1 ...
    0 0 0 -1 0 +1 +1 -1 +1 0 -1 -1 +1 0 0 +1 0 0 0 +1 0 0 0 0 +1 +1 0];                   % 24
    
    
    %% Length 91 - Table 15-7a
    [-1 0 +1 +1 +1 +1 -1 -1 +1 -1 -1 +1 -1 +1 +1 +1 +1 -1 +1 -1 -1 -1 +1 +1 -1 -1 +1 +1 ...
    +1 +1 +1 +1 -1 +1 +1 -1 +1 0 0 +1 -1 -1 +1 0 -1 -1 +1 0 +1 +1 +1 +1 +1 -1 -1 +1 ...
    +1 +1 -1 -1 0 -1 -1 0 +1 -1 +1 -1 -1 -1 -1 0 -1 +1 -1 +1 -1 +1 0 +1 -1 -1 +1 +1 -1 +1 ...
    -1 +1 +1 +1 0]; % 25
    
    [+1 +1 0 +1 -1 +1 -1 -1 -1 +1 +1 +1 +1 +1 -1 +1 -1 +1 +1 -1 -1 +1 -1 -1 +1 +1 -1 -1 -1 ...
    +1 -1 0 +1 +1 +1 0 -1 +1 +1 +1 +1 -1 +1 0 +1 0 -1 -1 0 +1 -1 +1 +1 -1 +1 +1 +1 +1 +1 +1 ...
    -1 -1 +1 -1 +1 +1 0 0 +1 +1 +1 -1 -1 0 +1 -1 -1 -1 -1 -1 +1 -1 0 +1 -1 +1 -1 +1 -1 -1 -1]; % 26
    
    [+1 +1 +1 -1 -1 +1 +1 +1 -1 -1 -1 +1 -1 +1 -1 0 -1 +1 -1 -1 0 +1 +1 -1 +1 -1 +1 0 -1 +1 ...
    +1 +1 +1 +1 +1 +1 +1 +1 +1 -1 -1 +1 -1 -1 +1 +1 -1 +1 +1 0 +1 +1 -1 +1 -1 +1 -1 -1 +1  ...
    -1 -1 +1 +1 +1 -1 -1 -1 0 -1 +1 +1 +1 -1 0 +1 0 0 -1 -1 -1 +1 +1 -1 +1 -1 -1 0 -1 +1 +1 0]; % 27
    
    [+1 +1 +1 +1 +1 -1 -1 +1 +1 +1 -1 +1 +1 -1 -1 -1 +1 -1 +1 -1 -1 0 +1 +1 -1 -1 -1 +1 -1 +1 0 ...
    +1 -1 -1 -1 -1 -1 +1 0 +1 +1 +1 -1 -1 +1 -1 +1 -1 -1 +1 -1 +1 +1 -1 +1 +1 +1 +1 0 -1 0 -1 +1 ...
    +1 0 0 +1 -1 +1 +1 +1 -1 +1 +1 -1 +1 0 -1 +1 0 -1 -1 +1 -1 -1 -1 +1 +1 +1 0 +1]; % 28
    
    [+1 -1 0 -1 -1 +1 -1 +1 +1 -1 -1 0 +1 +1 0 0 +1 +1 -1 +1 +1 -1 -1 -1 -1 -1 +1 +1 +1 +1 ...
    +1 +1 -1 0 +1 -1 -1 +1 -1 +1 +1 -1 -1 +1 -1 +1 +1 +1 -1 -1 +1 +1 +1 +1 +1 +1 +1 -1 +1 +1 ...
    +1 0 +1 -1 +1 -1 0 -1 0 -1 +1 +1 -1 -1 -1 +1 0 -1 -1 -1 +1 +1 0 -1 +1 -1 +1 -1 +1 +1 -1]; % 29
    
    [-1 +1 +1 0 -1 -1 0 +1 +1 -1 0 0 -1 -1 +1 +1 -1 +1 +1 -1 +1 -1 -1 +1 +1 +1 +1 +1 -1 -1 ...
    -1 +1 +1 +1 -1 +1 -1 0 +1 -1 +1 -1 +1 0 +1 +1 +1 +1 +1 -1 +1 +1 +1 -1 +1 +1 +1 -1 +1 0 ...
    -1 -1 -1 -1 -1 -1 +1 +1 -1 +1 +1 -1 +1 0 -1 -1 -1 -1 +1 -1 +1 -1 0 +1 0 +1 -1 +1 +1 +1 -1]; % 30
    
    [-1 +1 -1 +1 +1 0 +1 +1 +1 -1 -1 +1 +1 -1 0 +1 +1 0 0 -1 -1 +1 +1 -1 -1 +1 +1 +1 -1 +1 ...
    -1 -1 -1 -1 -1 -1 0 +1 +1 +1 -1 +1 +1 +1 +1 +1 -1 -1 +1 -1 +1 +1 -1 -1 -1 +1 -1 +1 -1 ...
    -1 -1 +1 -1 +1 0 -1 +1 -1 -1 0 +1 0 +1 -1 +1 +1 -1 +1 +1 0 +1 -1 +1 -1 -1 0 +1 +1 +1 +1 +1]; % 31
    
    [-1 +1 +1 +1 +1 +1 +1 +1 +1 +1 +1 +1 -1 -1 -1 +1 -1 +1 +1 -1 -1 +1 +1 0 0 -1 +1 -1 +1 0 ...
    -1 +1 -1 0 -1 +1 +1 -1 -1 +1 +1 +1 -1 +1 +1 +1 0 -1 -1 0 +1 +1 -1 +1 -1 +1 -1 0 -1 -1 -1 ...
    +1 +1 -1 0 -1 -1 -1 -1 +1 +1 +1 +1 -1 +1 -1 0 +1 0 -1 +1 -1 +1 +1 -1 +1 +1 -1 -1 +1 -1]; % 32
    };
end

code = Codes{codeIndex};


