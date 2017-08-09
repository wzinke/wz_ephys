function [Bdiff, Bse, bootdiff, bootresp1, bootresp2] = SPK_boot_RespDiff(mat1, mat2, nBoot, smplsz, npar)
% SPK_boot_RespDiff - calculate the bootstrapped difference of 2 SDF
% matrices
%  
% DESCRIPTION 
% 	This function calculates interspike intervals for 
% 	
%  
% SYNTAX 
%
% [Bdiff, Bste, diffmat] = SPK_boot_RespDiff(mat1, mat2, nBoot, smplsz)
%
%   Input:
%       mat1/mat2   2D matrix containing trial wise SDF for two conditions to be 
%                   compared. Rows correspond to trial number and collumns to time bin.
%                  
%       nBoot       number of bootstrap replications
%
%       smplsz      subsampling bootstrap that uses only <smplsz> of trials
%                   as bootstrap sample for each conditon.     
%            
%
%
%   Output:
%       Bdiff     
%
%
% REFERENCES 
%
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
%
% $Created : 31-Jul-2014 by wolf zinke
% $Modified: 

%  ======================================================================== 
%% prepare input data 

% number of bootstrap repetitions
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1000;
end

% number of bootstrap repetitions
if(~exist('smplsz','var') )
   smplsz = [];
end

% use non-parametric bootstrap instead of the parametric one
if(~exist('npar','var') || isempty(npar))
    npar = 0;
end

bootresp1 = SPK_boot_matrix(mat1, nBoot, smplsz, [], npar);
bootresp2 = SPK_boot_matrix(mat2, nBoot, smplsz, [], npar);

bootdiff = bootresp1 - bootresp2;

if(npar==0)
    Bdiff = nanmean(bootdiff,1);
    Bse   =  nanstd(bootdiff,[],1);
else
    Bdiff = nanmedian(bootdiff,2);
    Bse   =   prctile(bootdiff, [5 95], 1);
end


