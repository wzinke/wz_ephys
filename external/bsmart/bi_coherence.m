function coherence = bi_coherence(directory,n,fs)
% BI_COHERENCE Compute the pair coherence from the Bivariate models
% 
% Usage:
%   coherence = bi_coherence(directory,n,fs);
% 
% Input(s):
%   directory   - Onewindow_Coefficient or Movingwindow_Coefficient directory
%   n           - Number of frequency bins
%   fs          - Sampling rate
% 
% Example:
%   coherence = bi_coherence('\bsmartroot\Movingwindow_Coefficient',100,200)
% 
% See also: paircoherence.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 18-Sep-2007 14:49:09$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

path  = directory;
fil   = fullfile(path,'AR_C_*');
file1 = dir(fil);
c     = length(file1);
channel2 = (1+sqrt(1+8*c))/2;

% need Matlab symbolic toolbox
% TODO: better not to use this toolbox
% s     = sprintf('x^2-x-%d',c);
% dd    = solve(s);   
% dd    = double(dd);
% if dd(1)>0
%     channel2 = dd(1);
% else
%     channel2 = dd(2);
% end

pairco =[];
spect = [];
k = 0;
l = 1;
for i = 1:(channel2-1)
    for j = (i+1):channel2
        ii = num2str(i);
        jj = num2str(j);
        filAR = ['AR_C_' ii 'and' jj];
        fileAR = fullfile(path,filAR);
        filAN = ['AR_N_' ii 'and' jj];
        fileAN = fullfile(path,filAN);
        aredat = load(fileAR);
        % make it compatable with Matlab R14
        if size(aredat,2) == 1
            aredat = aredat.';    % if aredat is a column, change it to row
        end%if

        arndat = load(fileAN);
        % make it compatable with Matlab R14
        if size(arndat,2) == 1
            arndat = arndat.';    % if arndat is a column, change it to row
        end%if
        
        csd4d = MAR_csd4dmem(aredat,arndat,n,fs);
        [paircoh, partcoh, mulcoh, autospect] = MAR_coh3D(csd4d, fs);
        k = k+1;
        pairco(:,:,k) = paircoh(:,:,1); %#ok<AGROW>
        if i == 1
            l = l+1;
            spect(:,:,i) = autospect(:,:,1); %#ok<AGROW>
            spect(:,:,l) = autospect(:,:,2); %#ok<AGROW>
        end
    end
end

coherence = pairco;

end%bi_coherence

% [EOF]