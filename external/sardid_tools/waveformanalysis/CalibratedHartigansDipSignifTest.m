function [dip, p_value, xlow, xup, boot_dip]=CalibratedHartigansDipSignifTest(xpdf,nboot)

% Function [dip,p_value,xlow,xup,boot_dip]=HartigansDipSignifTest(xpdf,nboot)

% calibrated Hartigan's Dip Test (direct conversion from GAUSS code hpr.pdf)
%
% calculates Hartigan's DIP statistic and its significance for the empirical p.d.f  XPDF (vector of sample values)
% This routine calls the matlab routine 'HartigansDipTest' that actually calculates the DIP
% NBOOT is the user-supplied sample size of boot-strap
% Code by F. Mechler (27 August 2002)

% % calculate the DIP statistic from the empirical pdf
% [dip, xlow, xup, ifault, gcm, lcm, mn, mj]=HartigansDipTest(xpdf);
% N=length(xpdf);
%
% % calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
% boot_dip=[];
% for i=1:nboot
%    unifpdfboot=sort(unifrnd(0,1,1,N));
%    [unif_dip]=HartigansDipTest(unifpdfboot);
%    boot_dip=[boot_dip; unif_dip];
% end;
% boot_dip=sort(boot_dip);
% p_value=sum(dip<boot_dip)/nboot;

n = length(xpdf);
xvar = sort(xpdf);

%g = number of kernel evaluation points

g = 512;

% Create gamma values for student t and beta distributions

start = 1;
jump1 = .0005;
jump2 = .1;

% Generate values for Symmetric Beta Distribution

betavb = zeros(3800,1);
gammvb = zeros(3800,1);

for i=1:3000
    betavb(i) = start+(i-1)*jump1;
    gammvb(i) = 2^(4*betavb(i)-1)*(betavb(i)-1)*(beta(betavb(i),betavb(i)))^2;
end

start2 = start+2999*jump1;

for i=1:800
    t = start2+i*jump2;
    betavb(i+3000) = t;
    gammvb(i+3000) = 2^(4*t-1)*(t-1)*(beta(t,t))^2;
end

start1 = .505+realmin;
jump1 = .0005;
jump2 = .1;

% Generate values for Rescaled t Distribution

betavt = zeros(3800,1);
gammvt = zeros(3800,1);

for i=1:3000
    betavt(i) = start1+(i-1)*jump1;
    gammvt(i) = 2*betavt(i)*(beta(betavt(i)-.5, .5))^2;
end

start3 = start+2999*jump1;

for i=1:800
    t = start3 +i*jump2;
    betavt(i+3000) = t;
    gammvt(i+3000) = 2*t*(beta(t-.5, .5))^2;
end

% Transform the data, 1 mean; 2 sum; 3 logarithmic
%
% t = 0;
%
% %--------------------End of User ID section of Code---------------
%
% if (t == 1)
%     xvar = xvar/mean(xvar,2);
% elseif (t == 2)
%     xvar = xvar/sum(xvar,2);
% elseif (t == 3);
%     xvar = ln(xvar);
% end

[dip,xlow,xup]=HartigansDipTest(xpdf);

% Creating a kernel estimate for the unknown distribution factor

grid = xvar(1):(xvar(n)-xvar(1))/(g-1):xvar(n);
h = 1.06*std(xvar)*n^(-1/5);
h2 = ((4/7)^(1/9))*std(xvar)*n^(-1/9);

fhat = [];
for i=1:size(grid,2)
    z = (grid(i) - xvar)/h;
    kelx = ones(n,1).*exp(-0.5*(z.^2))/sqrt(2*pi);
    kelx = sum(kelx,1)/(size(xvar,1)*h);
    fhat = [fhat;kelx];
end
[fhat0,findx] = max(fhat);
xmax = grid(findx);
z = (xvar - xmax)/h2;
f2hat0 = ones(n,1).*((z.^2-1).*exp(-0.5*(z.^2))/sqrt(2*pi));
f2hat0 = sum(f2hat0)/(n*(h2^3));

dhat = abs(f2hat0)/(fhat0^3);

ds = zeros(nboot,1);
p_value = 0;

% Determine Betahat

% FIXME: this loop is very slow in Octave, try to vectorize it
for i=1:nboot
    if (dhat <= gammvb(3800))
        [~,track] = min(abs(gammvb-dhat));
        betahat = betavb(track);
        repvar = betarnd(betahat,betahat,n,1);
        d=HartigansDipTest(repvar);
        if (d > dip)
            p_value = p_value + 1;
        end
    elseif (dhat >= gammvt(3800))
        [~,track] = min(abs(gammvt-dhat));
        betahat = betavt(track);
        df = 2 * betahat - 1;
        repvar = trnd(df,n,1);
        d=HartigansDipTest(repvar);
        if (d > dip)
            p_value = p_value + 1;
        end
    else
        betahat = 2*pi;
        repvar = randn(n,1);
        d=HartigansDipTest(repvar);
        if (d > dip)
            p_value = p_value + 1;
        end
    end
    ds(i) = d;
end
p_value = p_value/nboot;

walpha  = [.90, .95, .99];
walpha  = quantile(ds, walpha);

boot_dip = ds;

%  fprintf('\n\n>>>>>>\tDip Statistic:\t\t%.3g\n',dip);
%  fprintf('\n>>>>>>\tp value:\t\t%.3g\n',p_value);
%  fprintf('\nW_alpha level quantiles:\t%.3g\t%.3g\t%.3g\n',walpha(1),walpha(2),walpha(3));
%  fprintf('\nModal Interval\t\t\t%.3g\t%.3g\n',xlow,xup);
%  fprintf('\n# of Observations:\t\t%.3g\n',n);
%  fprintf('\nD-Hat:\t\t\t\t%.3g\n',dhat);
%  fprintf('\nBeta-Hat:\t\t\t%.3g\n\n',betahat);
