function output = transres_RSA_dist_euclidean(decoding_out,chancelevel,cfg,data)

% function output = transres_RSA_dist_euclidean(decoding_out,chancelevel,cfg,data)
% 
% Calculates the euclidean distance between all datapoints of the full 
% datamatrix.
%
% 2013 Martin H.

tmp = sum(data.*data,2);
output = {sqrt(bsxfun(@plus,tmp,tmp')-(2*data)*data')};
