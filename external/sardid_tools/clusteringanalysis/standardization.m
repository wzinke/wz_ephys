function data = standardization(data)
% data --- observations x dimensions,
% standardization: every column is standardized within [0,1]

   nrow = size(data,1);
   colmin = min(data);
   colmax = max(data);
   dmax = colmax-colmin;
   data = data - repmat(colmin,nrow,1);
   data = data./repmat(dmax,nrow,1);
