function pairmat = wz_combpairs(vec)
% get all pairs of elements in a vector without replication
%
% wolf zinke, Nov 2014

vec = sort(unique(vec));
num_el = length(vec);

pairmat = nan(sum(1:num_el-1),2);

cnt = 1;
for(i=1:num_el-1)
    c_len = num_el - i;
    
    pairmat(cnt:cnt+c_len-1,1) = repmat(vec(i), c_len, 1);
    pairmat(cnt:cnt+c_len-1,2) = vec(i+1:end);
      
    cnt = cnt+c_len;
end

