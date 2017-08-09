function x = minor(x, i, j);

%  calculate a minor of matrix, X, corresponding to the element X_ij
% this is done by omitting the ith row and jth column 
%

% 
% Hualou Liang, 11/25/98, FAU
% 

[row, col, freq] = size(x);

if i>row | j>col
  error('the number of row or column is too large in calculating minor');
end

x(i,:,:) = [];
x(:,j,:) = [];

