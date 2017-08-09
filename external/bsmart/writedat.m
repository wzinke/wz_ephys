function writedat(fname,dat)
% WRITEDAT Write data in format that steve required
%   fname: file name of output file
%   dat: T x chans x Trial

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 16:03:28$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Hualou Liang
% revised 03/31/99, FAU

[T, CHN, TRLS] = size(dat);

fid = fopen(fname,'wb');
for i = 1:TRLS,
    fwrite(fid,squeeze(dat(:,:,i)),'float');
end%for
fclose(fid);

end%writedat

% [EOF]
