function valid_errorate(labels, truelabels)
% computing error rates for every clusters if true labels are given

nrow = length(truelabels);
nk = max(truelabels);
R=ind2cluster(labels);
S=ones(nrow,1);
high=[];
for i = 1:nk
   high(i)=mean(R{i});
end
[low, high]=sort(high);
for i = 1:nk
   low=high(i);
   low=R{low};
   S(low)=i;
   vtype=S(low)-truelabels(low);
   vtype=nonzeros(vtype);
   vtype=length(vtype);
   Rerror=100*vtype/length(low);
   fprintf('\n Error rate of cluster %d : %4.2f %%',i, Rerror);
end

S=S-truelabels;
S=nonzeros(S);
S=length(S);
Rerror=100*S/nrow;
fprintf('\n Error rate for all the data: %4.2f %% \n',Rerror);
