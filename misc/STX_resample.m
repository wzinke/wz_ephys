function [out_data, num_vals, smpl_sz] = STX_resample(in_data,DIM)
% Resampling Data with replacement.
%
% wolf zinke, 10.09.2006


if(exist('DIM','var') == 0)
    DIM = 1;
elseif(isempty(DIM) == 1)
    DIM = 1;
end

in_data = squeeze(in_data);

if(ndims(in_data) > 2)
    error('Function only works with 2-D datasets!');
end

data_dim = size(in_data);

if(DIM == 2)
    num_vals = data_dim(2);
    smpl_sz  = data_dim(1);
else
    num_vals = data_dim(1);
    smpl_sz  = data_dim(2);
end

% Generate random positions for the data
rndpos = round(1 + (smpl_sz-1) * rand(data_dim));

if(DIM == 2)
    for(i=1:num_vals)
        out_data(:,i) = in_data(rndpos(:,i),i);
    end
else
    for(i=1:num_vals)
        out_data(i,:) = in_data(i,rndpos(i,:));
    end
end


