function omat = nanRbind(mat1, mat2)
% nanRbind - concatenate rows and pad with NaN 
%
% DESCRIPTION 
%   Concatenate by row the two input variables that should be either a
%   row vector or a 2D matrix. 
%   A column vector will be interpreted row wise.
% 
% SYNTAX 
%   omat = nanRbind(mat1, mat2)
%
%   Input:
%       <mat1, mat>  Two row vectors or 2D matrices to combine row wise  
%
% wolf zinke

sz1 = size(mat1,2);
sz2 = size(mat2,2);

if(~isempty(mat1) && ~isempty(mat1))
    if(sz1 < sz2)
        mat1(:, sz1+1:sz2) = NaN;
    elseif(sz1 > sz2)
        mat2(:, sz2+1:sz1) = NaN;
    end
end

omat = [mat1; mat2];

