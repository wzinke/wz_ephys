function [bootmat, Bmean, Bse] = SPK_boot_matrix(inmat, nBoot, smplsz, DIM, npar)
% SPK_boot_matrix - calculate the bootstrapped estimate of the mean for a 2D matrix
% matrices
%
% DESCRIPTION
% 	This function resamples the rows (DIM=1, default) or columns (DIM=2) of
% 	a 2D matrix and calculates the bootstrap estimate of the mean and of
% 	the standard error (or median and mad if npar=1).
%
%  The first bootstrap sample will correspond to the raw data, unless
%  smplsz is used and different from the input size.
%
% WARNING: Be carefull when using large matrices and large numbers for bootstrap
%          repetitions, this could cause memory issues with this approach.
%          In this case it might
%
%
% SYNTAX
%
% [Bmean, Bste, bootmat] = SPK_boot_matrix(inmat, nBoot, smplsz, DIM)
%
%   Input:
%
%       inmat1    2D matrix containing trial wise SDF for two conditions to be
%                 compared. Rows correspond to trial number and collumns to time bin.
%
%       nBoot     number of bootstrap replications
%
%       smplsz    subsampling bootstrap that uses only <smplsz> of trials
%                 as bootstrap sample for each conditon.
%
%       DIM       Dimension to be bootstrapped. If DIM=1 (default)the rows will be
%                 resampled, if DIM=2, columns will be resampled
%
%       npar      use non-parametric bootstrap and determine median and MAD instead
%                 of mean and std.
%
%
%
%   Output:
%
%       Bmean
%
%       Bse
%
%       bootmat
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
%% check input and set default values

% number of bootstrap repetitions
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1000;
end

% what dimension of the matrix will be resampled
if(~exist('DIM','var') || isempty(DIM))
    DIM = 1;
end

% number of bootstrap repetitions
if(~exist('smplsz','var') || isempty(smplsz))
    smplsz = size(inmat, DIM);
end

% use non-parametric bootstrap instead of the parametric one
if(~exist('npar','var') || isempty(npar))
    npar = 0;
end


%  ========================================================================
%% do the calculation

% for simplicity of the code just transpose the matrix
if(DIM==2)
    inmat = t(inmat);
end

% determine the number of samples
N = size(inmat, 1);

% get a random matrix of sample indices
smplmat = randi(N, [smplsz, nBoot]);

% transform it to a third dimension, but make sure that rows are keeping
% the same index (should use reshape here, to make the matrix [Nrow, Ncol, Nboot],
% right now it is [Nrow, Nboot, Ncol]
% smplmat = repmat(smplmat, 1, 1, size(inmat,2));

% get the bootstrap samples
% allsmpl = inmat(smplmat);
bootmat = nan(nBoot, size(inmat, 2));

parfor(b=1:nBoot)
    csmpl = inmat(smplmat(:,b),:);
    if(npar==0)
       bootmat(b,:) = nanmean(csmpl);
    else
       bootmat(b,:) = nanmedian(csmpl);
    end
end

% calculate the bootstrap statistics
if(npar==0)
%     bootmat = squeeze(nanmean(allsmpl,1));
    if(nargout>1)
        Bmean = nanmean(bootmat,1);
    end
    if(nargout>2)
        Bse = nanstd(bootmat,[],1);
    end
else
%     bootmat = squeeze(nanmedian(allsmpl,1));
    if(nargout>1)
        Bmean = nanmedian(bootmat,2);
    end
    if(nargout>2)
        Bse = prctile(bootmat, [5 95], 1);
    end
end

% clear allsmpl

% if the matrix was transposed, do it again
if(DIM==2)
    Bmean   = t(Bmean);
    Bse     = t(Bse);
    bootmat = t(bootmat);
end

