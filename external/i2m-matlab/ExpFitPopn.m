BinSize = 15;
func = @(x,t)exp(-t./x);
%idx2 = find(lag2 >= 0 & lag2 < 20);
%[x y]=lsqcurvefit(func,500,lag2(idx2)*BinSize,tabPopn2(idx2),0,1000)
idx2 = find(lag2 >= 0 & lag2 < 100);
[x y]=lsqcurvefit(func,500,lag2(idx2)*BinSize,tabPopn2(idx2),0,1000)



plot(lag2(idx2)*BinSize,tabPopn2(idx2))
hold on
plot(lag2(idx2)*BinSize,func(x,lag2(idx2)*BinSize),'r')