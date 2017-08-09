function [resid] = whiteness_test(dat,window,order)
% WHITENESS_TEST Whiteness test
%
%   dat: data set in Matlab format
%   window: window length
%   order: model order
%
%   resid: residuals probabilities
%
% Example:
%   [resid] = whiteness_test(data,10,5)

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 22:43:55$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% parse directory
cdir = pwd;                 % find current directory
p = mfilename('fullpath');
fdir = fileparts(p);        % find function directory
cd(fdir);                   % change current dir to function dir

[points channel trail] = size(dat);
save channel channel -ascii;
save trail trail -ascii;
save points points -ascii;
save window window -ascii;
save order order -ascii;

writedat('dataset.bin',dat);
%!make -f Makefile.whiteness
if ispc
    status = eval(['unix ' '(''' 'opsswhite ' 'dataset.bin ' ' A ' 'Ve' ''')']);
else
    status = eval(['unix ' '(''' './opsswhite ' 'dataset.bin ' ' A ' 'Ve' ''')']);
end%if
% clean job one
delete channel
delete trail
delete points
delete window
delete order
delete dataset.bin
% if error
if status ~= 0
    error('Cannot test whiteness with ''opsswhite''!');
end%if
% otherwise
resid = load('resid.out');
% clean job two
delete resid.out
delete A
delete Ve

% restore directory
cd(cdir);

end%whiteness_test

% [EOF]