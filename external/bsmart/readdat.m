function dat = readdat(fname,NGPTS,nchns,ngdtrls)
% READDAT Read data in binary format and convert them into Matlab data
%
% Syntax:
%   dat = readdat(fname,NGPTS,nchns,ngdtrls);
%
% Input(s):
%   fname   -   Binary data file path and name
%   NGPTS   -   Number of points to be read in each trial
%   nchns   -   Number of channels to be read
%   ngdtrls -   Number of trials to be read
%
% Output(s):
%   dat     -   Three dimensional Matlab data array
%               = Points x Channels x Trials
%
% Example:
%   dat = readdat('test71.bin',18,15,137);
%
% See also: writedat.

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 16:03:28$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang
% revised 03/31/99, FAU

fid = fopen(fname,'rb','ieee-le');
dat = fread(fid,NGPTS*nchns*ngdtrls,'float');
dat = reshape(dat,[NGPTS nchns ngdtrls]);
fclose(fid);

end%function

% [EOF]