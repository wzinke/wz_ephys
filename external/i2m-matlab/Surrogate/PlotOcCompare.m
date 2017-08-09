function Djs=PlotOcCompare(Pb1,Pb2,Featureb,BinSize,TempSize,N,N_r,n_sample)

NormalisRate=1/(0.001*BinSize);
PointColor{1}='b.';
PointColor{2}='r.';
PointColor{3}='g.';
PointColor{4}='y.';
PointColor{5}='m.';
PointColor{6}='c.';

idx = find(Pb1>0 & Pb2 >0);
P1 = Pb1(idx);
P2 = Pb2(idx);
Feature = Featureb(idx);

MaxFeature=max(Feature)
% figure;

c=colormap;

for i=1:MaxFeature
    indF=find(Feature==i);
    if i==2
        hold on
    end
    loglog(NormalisRate*P1(indF),NormalisRate*P2(indF),'.','Color',c(ceil(size(c,1)*i/MaxFeature),:),'MarkerSize',20)
    hold on
end
a=[min(NormalisRate*P1) max(NormalisRate*P1)];
b=logspace(log10(min(P1)),log10(max(P1)),500);
y1 = NormalisRate*b + NormalisRate*1.96.*sqrt(b.*(1-b)./(N_r-1));
y2 = NormalisRate*b - NormalisRate*1.96.*sqrt(b.*(1-b)./(N_r-1));
loglog(NormalisRate*b,y1,'k--','Linewidth',2)
loglog(NormalisRate*b,y2,'k--','Linewidth',2)
loglog(a,a,'k','Linewidth',2);
% colorbar()
caxis([0 MaxFeature])
set(gca,'FontSize',14);
xlabel('Original stationnary rate (/s)','FontSize',16);
ylabel('Generated stationnary rate (/s)','FontSize',16);

Djs=BootStrapDJS(P1,P2,n_sample,N_r);

if nargin >4
   title(['Size of the patterns: ' int2str(N) ' x ' int2str(TempSize) '. DJS=' num2str(Djs(1))],'FontSize',16);
end
