clear, close,clc,

% make handles global so all function will 'see' all parameters they need
global handles 

% below set the requiered parameters as you wish
handles.binsize = 1; % window length used to bin data
handles.nscales = 4; % number of scales for wavelet decomposition
handles.nsurr = 50; % number of surrogates for computing shuffling distribution
handles.percentile = 95; % percentile of surrogate distribution for significance
handles.maxwvcoefs = 25; % maximum number of coefs to use
handles.minwvcoefs = 2; % mininum number of coefs to use

% below a mat file andd add the data to the global variable
load('simdata_ex2.mat')
handles.spiketimes = spiketimes; 
handles.class_id = class_id; 

WIfunc_binmatrix() % constructs binned matrix
WIfunc_wavedec() % computes wavelet decomposition
WIfunc_decode_5050() % classifies with 50/50 crossvalidation
% WIfunc_decode_leaveoneout() % classifies with leave-one out crossvalidation

%%

handles.matrices
handles.decode