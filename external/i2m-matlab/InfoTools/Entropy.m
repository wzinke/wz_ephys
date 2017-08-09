function S=Entropy(P)
%Brute force estimation
%P1=P/sum(P);%not sure about that...
S=-sum(P.*log(P));
