% Kaijun WANG, sunice9@yahoo.com, Oct. 2006, March 2007
function SR=validity_Index_mod(data,classlabel,N1,N,truelabels,nk,N2,Rd,Dist,dmax,nrow)


NC = N1:N;
labels = classlabel;
Rand = zeros(1,N);
Mirkin = zeros(1,N);
Hubert = zeros(1,N);
Sil = zeros(1,N);
DB = zeros(1,N);
CH = zeros(1,N);
Ha = zeros(1,N);
Hom = zeros(1,N);
Sep = zeros(1,N);

% (1) External validity indices when true labels are known
for i = NC
  [~,Rand(i), Mirkin(i), Hubert(i)] = ...
      valid_RandIndex(labels(:,i),truelabels);
end
if nk > 1
  valid_errorate(labels(:,nk), truelabels);   % error rate if true labels are given
end

Re = strcmp(Rd, 'euclidean');
% (2) Internal validity indices when true labels are unknown
for i = NC
   R = silhouette(data, labels(:,i), Rd);
   Sil(i) = mean(R);        % average Silhouette
   % Davies-Bouldin, Calinski-Harabasz
   [DB(i), CH(i), ~, Ha(i), ST] = valid_internal_deviation(data,labels(:,i), Re);
   S = ind2cluster(labels(:,i));
   [Hom(i), Sep(i)] = valid_internal_intra(Dist, S, Re, dmax);
end

%  Homogeneity-Separation
Hom = dmax*Hom;
Sep = dmax*Sep;

SR = [Rand; Mirkin; Hubert; Sil; DB; CH; Ha; Hom; Sep];
