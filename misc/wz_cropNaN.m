function [out_mat, valcnt] = wz_cropNaN(in_mat, minval, maxval, tonan)
% wz_cropNaN - removes collumns cosisting of NaN only froma 2D matrix
%  
% DESCRIPTION 
% Remove values smaller than <minval> (inclusive) or larger than <maxval> (exclusive),
% left align all values, and crop columns that contain only NaN values.
%  
% SYNTAX 
% [out_mat, valcnt] = wz_cropNaN(in_mat, minval, maxval, tonan)

% INPUT
%       in_mat  2D matrix to be cropped.
%
%       minval  Lower threshold of valid values. Only data that is larger
%               than <minval> will be included in the output matrix
%
%       maxval  Upper threshold of valid values. Only data that is smaller 
%               or equal than <minval> will be included in the output matrix
%
%       tonan   Value that will be interpreted as NaN prior cropping.

%
%   Output:
%       out_mat   cropped and aligned
%
%       valcnt    coeeficient of variation
%
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 31-Jul-2014 by wolf zinke
% $Modified: 23-Aug-2014 by wolf zinke: removed bug that reshaped values incorrectly


% ========================================================================= 
%% if specified set values outside the defined range to NaN
if(exist('minval','var') && ~isempty(minval))
    in_mat(in_mat < minval) = NaN;
end

if(exist('maxval','var') && ~isempty(maxval))
    in_mat(in_mat > maxval) = NaN;
end

if(exist('tonan','var') && ~isempty(tonan))
    in_mat(in_mat == tonan) = NaN;
end

% ========================================================================= 
%% remove columns that consists entirely of NaNs
vmat   = isfinite(in_mat);
valcnt = sum(vmat,2);

% ========================================================================= 
%% remove NaNs at the start of each row and correct for interspersed NaNs
out_mat = nan(size(in_mat,1), max(valcnt));

for(i=1:length(valcnt))
    out_mat(i,1:valcnt(i)) = in_mat(i,vmat(i,:));
end

% ========================================================================= 
%% below the snipet is buggy, apparently not reshaping the value vector correctly
% pmat = repmat(1:max(valcnt),size(in_mat,1),1);
% pmat = bsxfun(@le, pmat, valcnt(:));
% 
% out_mat(pmat) = in_mat(vmat);
