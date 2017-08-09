function [A,Ve] = one_mul_model(dat,order,startp,window)
% ONE_MUL_MODEL One window for multivariate model
% 
% Usage:
%   [A,Ve] = mov_mul_model (dat,order,startp,window);
% 
% Input(s):
%   dat: input file in matlab format;
%   order: model order
%   startp: start position of the window
%   window: window length
% 
% Output(s):
%   A   -   AR coefficient file
%	Ve  -   AR noise file

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 22:58:31$
% SHIS UT-Houston, Houston, TX 7730, USA.
% 
% Lei Xu, Hualou Liang

% parse directory
cdir = pwd;                 % find current directory
p = mfilename('fullpath');
fdir = fileparts(p);        % find function directory
cd(fdir);                   % change current dir to function dir

% processing
channel = size(dat,2);
trail = size(dat,3);
save channel channel -ascii;
save trail trail -ascii;
save points window -ascii;
save order order -ascii;

start = startp;
dat = dat(start:start+window-1,:,:);
writedat('dataset.bin',dat);

if ispc
    eval(['unix ' '(''' 'opssfull ' 'dataset.bin ' ' A ' 'Ve ' 'AIC' ''')'])
else
    eval(['unix ' '(''' './opssfull ' 'dataset.bin ' ' A ' 'Ve ' 'AIC' ''')']);
end%if
A = load('A');
Ve = load('Ve');

delete A
delete Ve
delete channel
delete trail
delete points
delete AIC
delete order
delete dataset.bin

% restore directory
cd(cdir);

end%one_mul_model

% [EOF]