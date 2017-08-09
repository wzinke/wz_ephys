function index=IndexPattern(rexp,rth,PexpTab,PthTab,BinSize)

NormalisRate=1/(0.001*BinSize);

a=(NormalisRate*PexpTab-rexp).^2 + (NormalisRate*PthTab-rth).^2;
ind=find(a==min(a));
index=ind(1);
