function dist=DJS(P,Q)
%Brute force estimation
M=0.5*(P+Q);

dist=0.5*(KL(P,M)+KL(Q,M));
