function omat = nanRbind(mat1, mat2)
% concatenate rows of two 2D matrices and fill with NaN if needed.
%
% wolf zinke

sz1 = size(mat1,2);
sz2 = size(mat2,2);

if(sz1 < sz2)
    mat1(:, sz1+1:sz2) = NaN;
elseif(sz1 > sz2)
    mat2(:, sz2+1:sz1) = NaN;
end

omat = [mat1; mat2];

