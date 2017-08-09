function dist=KL(P,Q)
%Brute force estimation
dist=sum(P.*log(P./Q));
