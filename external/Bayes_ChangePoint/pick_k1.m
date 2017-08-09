function num_chgpts = pick_k1(k)
% Function will sample from any probability distribution, initially created
% to sample a number of change points

dist = 0;
for i=1:length(k);  
    k(i)=dist+k(i); % Cumulative Probabilities (CDF)
    dist=k(i);
end

r=rand();           % Uniform random number between 0 and 1
i=1;
while (r>=k(i) && i<=length(k))
    i=i+1;          % Sample from cumulative probabilities (CDF)
end
num_chgpts=i;     
       