function pattern=ShowPattern(rexp,rth,PexpTab,PthTab,PatternList,BinSize)
index=IndexPattern(rexp,rth,PexpTab,PthTab,BinSize);
pattern=(PatternList(:,:,index)+1)/2;
