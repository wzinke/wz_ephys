function Djs=PlotOcFit(Pexp,Pth,Feature,BinSize,TempSize,N,N_r)
%Plot the Fig. 1 typical plot and estimate the Djs
%Pexp and Pth are the prediction and empirical estimation to compare
NormalisRate=1/(0.001*BinSize);
PointColor{1}='b.';
PointColor{2}='r.';
PointColor{3}='g.';
PointColor{4}='y.';
PointColor{5}='m.';
PointColor{6}='c.';

MaxFeature=max(Feature)
figure;

c=colormap;

for i=1:MaxFeature
    indF=find(Feature==i);
    if i==2
        hold on
    end
    loglog(NormalisRate*Pexp(indF),NormalisRate*Pth(indF),'.','Color',c(ceil(size(c,1)*i/MaxFeature),:),'MarkerSize',20)
    hold on
end
a=[min(NormalisRate*Pexp) max(NormalisRate*Pexp)];
b=logspace(log10(min(Pexp)),log10(max(Pexp)),500);
y1 = NormalisRate*b + NormalisRate*1.96.*sqrt(b.*(1-b)./(N_r-1));
y2 = NormalisRate*b - NormalisRate*1.96.*sqrt(b.*(1-b)./(N_r-1));
loglog(NormalisRate*b,y1,'k--','Linewidth',2)
loglog(NormalisRate*b,y2,'k--','Linewidth',2)
loglog(a,a,'k','Linewidth',2);
colorbar()
caxis([0 MaxFeature])
set(gca,'FontSize',14);
xlabel('Observed pattern rate (/s)','FontSize',16);
ylabel('Predicted pattern rate (/s)','FontSize',16);

Djs=DJS(Pexp,Pth);

if nargin >4
    title(['Size of the patterns: ' int2str(N) ' x ' int2str(TempSize) '. DJS=' num2str(Djs(1))],'FontSize',16);
end
