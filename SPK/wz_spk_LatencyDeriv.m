function [lat,fs,tfs] = wz_spk_LatencyDeriv(spk, timwin, bw)
% adapted from the estimateLatencyByDerivative function  created by
% Martin Nawrot (2207) and modified by  Ralph Meier for the FIND toolbox (2008).
% 
% References: 
% Meier et al. 2008, Neural Networks, Spec. Issue Neuroinformatics  


% Default Parameters:
m=1;

if(~exist('timwin','var') || isempty(timwin))
    timwin = spk.timewindow;
end

if(~exist('timwin','var') || isempty(timwin))
    timwin = spk.timewindow;
end

nl = -n;      % number of data points to the left
nr = -nl;     % number of data points to the rigth
ld = 1;       % 0: smooting, 1: numerical derivative
savgol = makeSavgol('n',n,'ld',ld,'m',m, ...
    'wf','Welch','h',1e-3);
le=length(savgol);

if Tdefault
    T=[1 length(s)-le];
end

fs = filter(fliplr(savgol),1,s);
%fs = filter((savgol),1,s);
fs=fs(le:end);
tfs=(1:length(fs))+nr;
[m,lat] = max(fs(T(1)+1:T(2)));
lat=lat(1)+T(1);
lat=lat+nr;