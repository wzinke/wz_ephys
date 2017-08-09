function power= bi_power(directory,n,fs)
% BI_POWER Compute the auto power from the Bivariate models
%
% Usage:
%   power = bi_power(directory,n,fs)
%
%   directory: Onewindow_Coefficient or Movingwindow_Coefficient directory
%   n: Number of frequency bins
%   fs: Sampling rate
%
% Example:
%   power = bi_power('bsmartroot\Movingwindow_Coefficient',100,200)

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.5$ $Date: 18-Sep-2007 19:11:23$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

path = directory;
fil = fullfile(path,'AR_C_*');
file1 = dir(fil);
c = length(file1);
channel2 = (1+sqrt(1+8*c))/2;

pairco=[];
spect=[];
k=0;l=1;
for i=1:(channel2-1)
    for j=(i+1):channel2
        ii=num2str(i);jj=num2str(j);
        filAR=['AR_C_' ii 'and' jj];
        fileAR=fullfile(path,filAR);
        filAN=['AR_N_' ii 'and' jj];
        fileAN=fullfile(path,filAN);
        aredat=load(fileAR);
        % make it compatable with Matlab R14
        if size(aredat,2) == 1
            aredat = aredat.';    % if aredat is a column, change it to row
        end%if
        arndat=load(fileAN);
        % make it compatable with Matlab R14
        if size(arndat,2) == 1
            arndat = arndat.';    % if arndat is a column, change it to row       end%if
        end
        csd4d = MAR_csd4dmem(aredat,arndat,n,fs);
        [paircoh, partcoh, mulcoh, autospect] = MAR_coh3D(csd4d, fs);
        k=k+1;
        pairco(:,:,k)=paircoh(:,:,1); %#ok<AGROW>

        if i==1
            l=l+1;
            spect(:,:,i)=autospect(:,:,1); %#ok<AGROW>
            spect(:,:,l)=autospect(:,:,2); %#ok<AGROW>
        end
    end
end

% change into power
power = abs(spect).^2;

end%bi_power

% [EOF]