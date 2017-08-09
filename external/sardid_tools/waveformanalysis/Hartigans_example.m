% hartigans dip test demo

% create some obviously unimodal and bimodal Gaussian distributions just to
% see what dip statistic does
% Nic Price 2006

clear;
clc;
close all;

Center = 0.1:0.025:0.3;
sig = 0.1;
nboot = 10000;

npoints1 = 3*80;
npoints2 = 3*320;

xbin=-0.8:0.025:0.8;
for a = 3:5 % 1:length(Center)
  xpdf1 = -Center(a)+sig*randn(1,npoints1);
  xpdf2 = Center(a)+sig*randn(1,npoints2);
  xpdf(a,:) = sort([xpdf1 xpdf2]);
  [dip(a), p(a)] = HartigansDipSignifTest(xpdf(a,:), nboot);
  [dipCal(a), pCal(a)] = CalibratedHartigansDipSignifTest(xpdf(a,:)', nboot);

  figure('color','none','visible','off')
  hold on
  set(gca,'layer','top','color','none')
  n=hist(xpdf(a,:),xbin);
  h=bar(xbin,n,1,'LineWidth',2);
  set(h,'FaceColor',[0.8 0.8 0.8],'BarWidth',1);
  % plot([median(xpdf1) median(xpdf2)],[.9*max(n),1.1*max(n)],'kv','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2)
  plot([median(xpdf1) median(xpdf2)],[110,130],'kv','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2)
  % plot([mean(xpdf1)-1.96*sig mean(xpdf1)+1.96*sig;mean(xpdf2)-1.96*sig mean(xpdf2)+1.96*sig]',[.9*max(n)*ones(2,1),1.1*max(n)*ones(2,1)],'k-','LineWidth',2)
  % plot([quantile(xpdf1,[0.025 0.975]);quantile(xpdf2,[0.025 0.975])]',[.9*max(n)*ones(2,1),1.1*max(n)*ones(2,1)],'k-','LineWidth',2)
  plot([quantile(xpdf1,[0.025 0.975]);quantile(xpdf2,[0.025 0.975])]',[110*ones(2,1),130*ones(2,1)],'k-','LineWidth',2)
  if pCal(a) < 0.00095 && p(a) < 0.00095
    title(['H. p-val < 0.001, Hcal. p-val < 0.001'],'fontsize',24)
  elseif pCal(a) < 0.00095
    title(['H. p-val = ',num2str(p(a),3),', Hcal. p-val < 0.001'],'fontsize',24)
  else
    title(['H. p-val = ',num2str(p(a),3),', Hcal. p-val = ',num2str(pCal(a),3)],'fontsize',24)
  end
  % axis([-0.8 0.8 0 1.15*max(n)])
  axis([-0.6 0.6 0 140])
  xlabel('x','fontSize',24)
  ylabel('Count','fontSize',24)
  set(gca,'fontSize',24,'LineWidth',2,'TickDir','out','Box','off','color','none','Xtick',-.6:.3:.6,'Ytick',0:40:140)
  filename = ['HartigansExamplePanels/HartigansRandnExample_',num2str(a),'.svg'];
  plot2svg(filename)
end
