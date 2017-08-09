function [dip,cal_H_p_val_dip,xl,xu,ind_nar,ind_bro,ind_fuz,aic_1,aic_2,bic_1,bic_2]=waveformSeparation(X,Xname,Xlabel,Xlim,mu,sigma,pcomponents,FIGSDIR)
% function [dip,pdip,xl,xu,aic_1,aic_2,bic_1,bic_2]=waveformSeparation(X,Xname,Xlabel)
% given the ARGUMENTS:
% X: waveform measurements
% Xname: waveform measure name
% Xlabel: the xlabel in figures
% Xlim: limits of the x-axis
% mu: initial mean of each component in the 2-Gaussian fit
% sigma: initial covariance matrix of each component in the 2-Gaussian fit
% pcomponents: initial mixing proportions of each component in the 2-Gaussian fit
% this function:
% a) computes Calibrated Hartigans' dip test.
% OUTPUTS: [dip,pdip,xl,xu]
% b) plots its optimized histogram (depends on sshist)
% c) separates X into three groups: narrow, fuzzy and broad neurons according to a 2-Gaussian model
% OUTPUTS: [ind_nar,ind_bro,ind_fuz]
% d) plots colored histogram for the three groups
% e) computes Akaike's and Bayesian criteria of 1- and 2- Gaussian model to determine which is a better model, after taking into account the different number of parameters in the two fits.
% OUTPUTS: [aic_1,aic_2,bic_1,bic_2]

% Hartigan Dip Test for unimodality
nboot = 10000; % sample size of bootstrap
% FIXME: unnecessary loops in all Hartigans' related functions due to direct conversion from GAUSS code, this code needs to be vectorized
[dip,H_p_val_dip,xl,xu,boost_dip] = HartigansDipSignifTest(X, nboot);

% optimized histogram
xmin = min(boost_dip);
xmax = max(boost_dip);
range = xmax-xmin;
nbins = sshist(boost_dip);
step = range/nbins;
xbin = xmin-step:step:xmax+step;
n = hist(boost_dip,xbin);

% plotting: bootstrap sample of size nboot of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
hb=bar(xbin,n,'linewidth',1);
set(hb,'EdgeColor','k','FaceColor',[0.8 0.8 0.8],'BarWidth',1);
hp=plot([dip dip],[0 max(n)*1.05],'r','LineWidth',2);
xlabel('Dip','fontSize',16)
ylabel('Count','fontSize',16)
xlim([min([xmin,dip])-step,max([xmax,dip])]+step)
hl = legend([hb,hp],{'Bootstrap sample',['Dip ',Xlabel]},'Location','Best');
set(hl,'Visible', 'Off');
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
title(['dip = ',num2str(dip,2), ', p-val = ',num2str(H_p_val_dip,2)],'fontsize',16)
plot2svg([FIGSDIR,'/',Xname,'_OrigHart.svg']);

[dip,cal_H_p_val_dip,xl,xu,boost_dip] = CalibratedHartigansDipSignifTest(X, nboot);
% optimized histogram
xmin = min(boost_dip);
xmax = max(boost_dip);
range = xmax-xmin;
nbins = sshist(boost_dip);
step = range/nbins;
xbin = xmin-step:step:xmax+step;
n = hist(boost_dip,xbin);

% plotting: bootstrap sample of size nboot of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
hb=bar(xbin,n,'linewidth',1);
set(hb,'EdgeColor','k','FaceColor',[0.8 0.8 0.8],'BarWidth',1);
hp=plot([dip dip],[0 max(n)*1.05],'r','LineWidth',2);
xlabel('Dip','fontSize',16)
ylabel('Count','fontSize',16)
xlim([min([xmin,dip])-step,max([xmax,dip])]+step)
hl = legend([hb,hp],{'Bootstrap sample',['Dip ',Xlabel]},'Location','Best');
set(hl,'Visible', 'Off');
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
title(['dip = ',num2str(dip,2), ', p-val = ',num2str(cal_H_p_val_dip,2)],'fontsize',16)
plot2svg([FIGSDIR,'/',Xname,'_CalibHart.svg']);

% optimized histogram
xmin = min(X);
xmax = max(X);
range = xmax-xmin;
nbins = sshist(X);
step = range/nbins;
xbin = xmin-step:step:xmax+step;
n = hist(X,xbin);
area = trapz(xbin,n);

% Fits of 1- and 2- Gaussian model
% FIXME: gmdistribution.fit is not available in GNU Octave

% mu and sigma estimation
options = statset('MaxIter',1000,'TolFun',1e-12);
f_1 = gmdistribution.fit(X,1,'Options',options); % single component, unimodal

S.mu=mu;
S.Sigma=sigma;
S.Sigma=reshape(S.Sigma,1,1,2);
S.PComponents=pcomponents;
f_2 = gmdistribution.fit(X,2,'Options',options,'Start',S); % two components, bimodal

% Akaike's and Bayesian Information Criteria
% for a given criterion, the smaller the better
% Akaike
aic_1=f_1.AIC;
aic_2=f_2.AIC;
% Bayes
bic_1=f_1.BIC;
bic_2=f_2.BIC;

% generating the fits
npoints=1000;
xfit = xmin:(xmax-xmin)/(npoints-1):xmax;
f_1a = normpdf(xfit,f_1.mu,sqrt(f_1.Sigma));
gauss_1=area*f_1a;

[~,id1]=min(f_2.mu);
[~,id2]=max(f_2.mu);
f_2a = normpdf(xfit,f_2.mu(id1),sqrt(f_2.Sigma(id1)));
f_2a_cum = normcdf(xfit,f_2.mu(id1),sqrt(f_2.Sigma(id1)));
f_2b = normpdf(xfit,f_2.mu(id2),sqrt(f_2.Sigma(id2)));
f_2b_cum = normcdf(xfit,f_2.mu(id2),sqrt(f_2.Sigma(id2)));

gauss_2=area*(f_2a*f_2.PComponents(id1)+f_2b*f_2.PComponents(id2));
gauss_2a=area*f_2a*f_2.PComponents(id1);
gauss_2b=area*f_2b*f_2.PComponents(id2);

% plotting

% optimized range
iU = find(gauss_2b>gauss_2a);
yU = gauss_2b(iU(1))-gauss_2a(iU(1));
xU = xfit(iU(1));
iL = find(gauss_2a>gauss_2b);
yL = gauss_2b(iL(end))-gauss_2a(iL(end));
xL = xfit(iL(end));
m = (yU-yL)/(xU-xL);
xintersect1 = -yL/m+xL;
xintersect2 = -yU/m+xU;
xintersect = 0.5*(xintersect1+xintersect2);
x1 = xintersect-0.5*step:-step:xmin-step;
x2 = xintersect+0.5*step:step:xmax+step;
xbin = sort([x1,x2]);
n = hist(X,xbin);

% all the histogram in monochrome
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h=bar(xbin,n,'LineWidth',2);
set(h,'EdgeColor','k','FaceColor',[0.8 0.8 0.8],'BarWidth',1);
xlabel([Xlabel,' (ms)'],'fontSize',24)
ylabel('Number of neurons','fontSize',24)
if H_p_val_dip < 0.00095 && cal_H_p_val_dip < 0.00095
  title(['H. p-val < 0.0005, Hcal. p-val < 0.0005'],'fontsize',24)
elseif cal_H_p_val_dip < 0.00095
  title(['H. p-val = ',num2str(H_p_val_dip,3),', Hcal. p-val < 0.001'],'fontsize',24)
else
  title(['H. p-val = ',num2str(H_p_val_dip,3),', Hcal. p-val = ',num2str(cal_H_p_val_dip,3)],'fontsize',24)
end
set(gca,'fontSize',24,'LineWidth',2,'TickDir','out','Box','off','XTick',-.3:.1:.7,'YTick',0:20:80)
xlim(Xlim)
plot2svg([FIGSDIR,'/',Xname,'.svg']);

% cutoffs
area_2a=cumtrapz(xfit,gauss_2a);
area_2a=area_2a(end);
area_2b=cumtrapz(xfit,gauss_2b);
area_2b=area_2b(end);
thresh_a = xfit((area_2a-cumtrapz(xfit,gauss_2a))<=10*cumtrapz(xfit,gauss_2b));
thresh_a = thresh_a(1);
thresh_b = xfit(10*(area_2a-cumtrapz(xfit,gauss_2a))<=cumtrapz(xfit,gauss_2b));
thresh_b = thresh_b(1);

% separation
ind_nar = find(X<=thresh_a);
ind_bro = find(X>=thresh_b);
ind_fuz = find(X>thresh_a & X<thresh_b);
y1=hist(X(ind_nar),xbin);
y2=hist(X(ind_bro),xbin);
y3=hist(X(ind_fuz),xbin);

% plotting: histogram of different types in different colors
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h = bar(xbin, y1,1,'LineWidth',1);
set(h,'EdgeColor','k','FaceColor','r','BarWidth',1);
hpatch = get(h,'children');
set(hpatch,'FaceAlpha',0.5);
h = bar(xbin, y2,1,'LineWidth',1);
set(h,'EdgeColor','k','FaceColor','b','BarWidth',1);
hpatch = get(h,'children');
set(hpatch,'FaceAlpha',0.5);
h = bar(xbin, y3,1,'LineWidth',1);
set(h,'EdgeColor','k','FaceColor','k','BarWidth',1);
hpatch = get(h,'children');
set(hpatch,'FaceAlpha',0.5);
plot(xfit,gauss_2,'k','linewidth',2)
plot(xfit,gauss_2a,'r','linewidth',2)
plot(xfit,gauss_2b,'b','linewidth',2)
xlabel([Xlabel,' (ms)'],'fontSize',16)
ylabel('Number of neurons','fontSize',16)
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off','XTick',-.3:.1:.7,'YTick',0:20:80)
xlim(Xlim)
plot2svg([FIGSDIR,'/',Xname,'_sep.svg']);
