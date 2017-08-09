function clust=ind2cluster(ind)

k=max(ind);
for i=1:k
  t=find(ind==i);
  clust{i}=t;
end