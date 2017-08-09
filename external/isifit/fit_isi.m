function info = fit_isi(allspikes, timediv, plotit)
%************** this will run maximum likelihood/ EM fitting
%************** on the spike trains ************************
%*** allspikes:  full spike distribution at 0.5 ms bin resolution
%*** timdiv:     time discreteness of spike times, 2 default =>(1/2) ms
%*** plotit:     if it is desired to show plots of the results
%output: info - a structure containing various results
%******** info.redflag - if not zero, there was an error (too few spikes?)    
%******** info.isibin - ISI frequency distribution
%******** info.isicnt - raw counts of times ISI was test (short ISI tested more)
%******** info.meanisi - mean ISI value
%******** info.zall - model fit of gaussian refractory * exponential decay
%                zall.wpar1 - [tau refractory, sigma refractory, rate]
%                zall.wpar2 - [same parameters, if 2nd exponent (not used)]
%                zall.p     - [ratio of mixtures of exponents, fixed to 1]
%                zall.varexp - variance explained of fit in 1st 50ms
%                zall.pfit  - actual fit values for full ISI distribution
%******* info.all - same as zall above, but with model including two
%*******            fit exponentials, the first modeling burst firing
%*******            and the second with regular rate firing


%*****************************************
global maxit Memb MembR Isi Fisi1 Fisi2 AFisi AFisi1 AFisi2 T WREF TP;
   
   TP = timediv;
   %************* insert in data to pass to routine for computing
   %************* the isi distribution (designed to take multiple arrays,
   %************* in case we wanted to fit attended vs unattended seperate)
   spikearray = cell(1,1);
   spikearray{1} = allspikes;
   
   %****************************************************************
   disp(sprintf('Beginning computation of ISI distribution'));
   [isibin,isicnt,redflag] = comp_isifreq(spikearray);
   
   if (plotit==1)
     figure;
     subplot(3,3,1);
     Xt = (1/TP) * (1:size(isibin,2));
     plot(Xt,isibin(1,:),'k.'); hold on;
     xlabel('Time (ms)');
     ylabel('Frequency');
     title('Raw ISI Distribution');
     subplot(3,3,2);
     Xt = (1/TP) * (1:size(isibin,2));
     x = find( Xt < 50 );  % find only first 50 ms to plot
     XX = 1:size(x,2);
     plot(Xt(XX),isibin(1,XX),'k.'); hold on;
     xlabel('Time (ms)');
     ylabel('Frequency');
     title('Raw ISI Dist (Zoom In)');
   end
   
   %**************** redflag any unit with no spikes in ISI testing
   info.redflag = redflag;
   if (redflag >= 1)
       disp(sprintf('\n Zero spike unit, no ISI analysis'));
       info = [];
       info.redflag = redflag;
       return;
   end
   
   %******** compute mean isi ***********************************
   %******* and store ISI distribution to info structure ********
   meanisi = 0;
   for i = 1:size(isibin,2)
       meanisi = meanisi + isibin(1,i)*i;  %set in 0.5ms bins
   end
   meanisi = (meanisi/TP);
   disp(sprintf('Mean Isi is %8.5f',meanisi));
   info.isibin = isibin;
   info.isicnt = isicnt;
   info.meanisi = meanisi;
   
   %*************** stop analysis if no spikes ************
   if (info.meanisi == 0)
      info = [];
      info.redflag = 1;   % not any spikes to measure rate?
      disp(sprintf('\n Zero spike unit, no ISI analysis'));
      return; 
   end
 
   AllIsi = isibin(1,:);
   
   %************** parameters to use in optimizations
   crit = 0.0000001;
   maxit = 150;
   M = size(isibin,2);
   T = M;  %use full set of possible measured isi's
   X = 1:200;
   %********** gradient descent maximization ***************
   opts = optimset('fminsearch'); %get default values
   optimset(opts,'MaxFunEvals',(10000*6));  % default was 200*6
   
   
   %********* throw a blanket over it for initial search ***
   %******* avoid local minima by fitting from many starting points ***
   
   %***** First model fit includes a sigmoid refractory (taur,tsig)
   %***** and a exponential decaying poisson or rate (a gaussian with parameters
   %***** taur and tsig convolved with an exponential distribution of rate)
   
   MN = 20;
   startgrid = zeros(MN,MN,3);
   for i = 1:MN
       taur = 0.1 + (0.2 * (1.4^(i-1)));   % time of taur
       tsig = taur/3;   % fixed tightness
       for j = 1:MN
           rate = 1.4^(j-1);   % 1 up to 400 Hz
           startgrid(i,j,1) = taur;
           startgrid(i,j,2) = tsig;
           startgrid(i,j,3) = rate;
       end
   end    
   best_err = 10000000;
   Isi = AllIsi;
   Isi = Isi / sum(sum( Isi ));  % must sum to one
   Memb = 0.5 * ones(1,T);
   errsurf = zeros(MN,MN);
   disp(sprintf('\n Running initial parameter search'));
   for i = 1:MN
       disp(sprintf('Fast fun %d of %d %10.5f',i,MN,best_err));
       for j = 1:MN
            taur = startgrid(i,j,1);
            tsig = startgrid(i,j,2);
            rate = startgrid(i,j,3);
            wparm1 = [log(taur),log(tsig),log(rate)];
            Fisi1 = isifun(wparm1);
            Fisi2 = Fisi1;
            errb = isicost3(0);
            errsurf(i,j) = errb;
            if (errb < best_err)
                taurbest = taur;
                tsigbest = tsig;
                ratebest = rate;
                best_err = errb;
            end
       end
   end
   %************* if plotit > 1 then show error surface *********
   if (plotit>1)
       subplot(3,3,3);
       imagesc(errsurf,[min(min(errsurf)) max(max(errsurf))]);
       colorbar;
       xlabel('Tau of Refractory');
       ylabel('Sigma of Refractory');
       title('Initial Parameter Search');
   end
   
   %******** Now that the blanket search found a good starting point,
   %******** perform an actual gradient descent search from that point.
   %************ used best params of blanket search for gradient search ***
   info.zall = single_fit(AllIsi,AllIsi,[log(taurbest),log(tsigbest),log(ratebest)],...
                               plotit,1,1);
   %*************** these fits allow the refractory to change also

   
   %********** Now we will fit a second model, this one also has a 
   %********** refractory period defined by a gaussian, but then instead
   %********** of drawing from a single exponential decay rate, it draws
   %********** from a mixture of two exponential decays, with probability
   %********** p from one with a higher rate (representing burst events)
   %********** (1-p) from one with a regular rate 
   
   %*************** same strategy could be applied in the 2-dual fit?
   %********* throw a blanket over it for initial search ***
   %******* avoid local minima by fitting from many starting points ***
   MN = 8;
   startgrid = zeros(MN,MN,MN,MN,MN,7);
   for i = 1:MN
       taur = 0.1 + (0.2 * (1.4^(i+1)));   % time of 0.8 up to 3.2
       tsig = taur/3;   % fixed tightness
       for j = 1:MN
           rate = 3^j;   % 4 up to 1024 Hz
           for ii = 1:MN
               taura = 0.1 + (0.2 * (2.5^(ii)));   % 3.2 up to 110
               tsiga = taura/3;
               for jj = 1:MN
                   ratea = 2^(jj-1);  %1 up to 40
                   for pp = 1:MN
                       p = 0.1 + ((pp-1)*0.08);  %0.1 up to 0.9
                       startgrid(i,j,ii,jj,pp,1) = taur;
                       startgrid(i,j,ii,jj,pp,2) = tsig;
                       startgrid(i,j,ii,jj,pp,3) = rate;
                       startgrid(i,j,ii,jj,pp,4) = taur; %taura;
                       startgrid(i,j,ii,jj,pp,5) = tsig; %tsiga;
                       startgrid(i,j,ii,jj,pp,6) = ratea;
                       startgrid(i,j,ii,jj,pp,7) = p;
                   end %pp
               end %jj
           end %ii
       end %j
   end %i
   %***************** selected many many possible starting points, test
   %them now  to find the best one ************************************
   best_err = 10000000;
   Isi = AllIsi;
   Isi = Isi / sum(sum( Isi ));  % must sum to one
   Memb = 0.5 * ones(1,T);
   errsurf = zeros(MN,MN,MN,MN,MN);
   disp(sprintf('\n Running initial parameter search, 2nd model'));
   tic
   for i = 1:MN
       for j = 1:MN
           disp(sprintf('Fast fun %d %d of %d %10.5f',i,j,MN,best_err));
       for ii = 1:MN
       for jj = 1:MN
       for pp = 1:MN
            taur = startgrid(i,j,ii,jj,pp,1);
            tsig = startgrid(i,j,ii,jj,pp,2);
            rate = startgrid(i,j,ii,jj,pp,3);
            taura = startgrid(i,j,ii,jj,pp,4);
            tsiga = startgrid(i,j,ii,jj,pp,5);
            ratea = startgrid(i,j,ii,jj,pp,6);
            pref = startgrid(i,j,ii,jj,pp,7);
            wparm1 = [log(taur),log(tsig),log(rate)];
            wparm2 = [log(taura),log(tsiga),log(ratea)];
            Fisi1 = isifun(wparm1);
            Fisi2 = isifun(wparm2);  
            OP = log( ((1-pref)/pref) );
            errb = isicost3(OP);
            errsurf(i,j) = errb;
            if (errb < best_err)
                taurbest = taur;
                tsigbest = tsig;
                ratebest = rate;
                taurabest = taura;
                tsigabest = tsiga;
                rateabest = ratea;
                pbest = pref;
                best_err = errb;
            end
       end %pp
       end %jj
       end %ii
       end %j
   end
   
   %*************** do dual poisson fits, refractory free to change ***********
   info.all = dual_fit(AllIsi,AllIsi,[log(taurbest),log(tsigbest),log(ratebest)],...
                   [log(taurabest),log(tsigabest),log(rateabest)],pbest,plotit,1,0,1);
   if (info.all.wparm1(3)<info.all.wparm2(3))  % make sure high fire in first position
       tmp = info.all.wparm1;
       info.all.wparm1 = info.all.wparm2;
       info.all.wparm2 = tmp;
       info.all.p = 1 - info.all.p;
   end
   
return;



%*************************************************
function ito = single_fit(fitisi,tisi,awref,plotit,plotnum,withref)
%******** fit recovery function and poisson isi dist
 global maxit Memb MembR Isi Fisi1 Fisi2 AFisi AFisi1 AFisi2 T WREF TP;

   %******************************************
   crit = 0.0000001;
   opts = optimset('fminsearch'); %get default values
   optimset(opts,'MaxFunEvals',(10000*6));  % default was 200*6
   
   %*********************
   disp(sprintf('\n Running gradient descent from best guess'));
    
   %***************************************
   XX = 1:size(fitisi,2);
   Xt = (1/TP) * XX;
   x = find( Xt < 50 );  % find only first 50 ms to plot
   X = 1:size(x,2);
   
   %*******************************************************
   dist = 10000000;
   wref = [awref(1),awref(2)];
   wparm1 = [wref,awref(3)];    % single time scale
   wparm2 = wparm1;    % don't fit longer time scale
   Fisi1 = isifun(wparm1);
   Fisi2 = isifun(wparm2);
   
   Isi = fitisi(1:T);
   Isi = Isi / sum(sum(Isi));
   OP =  0;
   P = 0.5; % fixed to 0.5, and tau1 is fixed to equal tau2
   [Memb,MembR] = membifun(wparm1,wparm2,P);  % initial comp of member func
   oerr = isicost3(OP);  % initial estimate of error
   iter = 0;
   while (dist>crit) & (iter<maxit)
       
      %********** fit wparm1 ***************
      WREF = [wref,wparm1(3),wparm2(3),1];
      wpar1 = wparm1(3);
      [wpar1,err1] = fminsearch( @isicost1a, wpar1, opts);
      WREF(3) = wpar1;
      wpar2 = wpar1;
      WREF(4) = wpar2;
      if (withref == 1)
        [wref,err3] = fminsearch( @isicost1c, wref, opts);
        WREF(1) = wref(1);
        WREF(2) = wref(2);
      end
      wparm1(1) = wref(1);
      wparm1(2) = wref(2);
      wparm1(3) = wpar1;
      wparm2(1) = wref(1);
      wparm2(2) = wref(2);
      wparm2(3) = wpar2;
      %*************************************
      Fisi1 = isifun(wparm1);
      Fisi2 = isifun(wparm2);
      %********** fit membership ***********
      errb = isicost3(OP);
      P = 0.5; %1/(1+exp(OP));
      [Memb,MembR] = membifun(wparm1,wparm2,P);
      %*************************************
      dist = (abs(oerr - errb)/abs(oerr)) * 100;
      iter = iter + 1;
      oerr = errb;
      %************************************
      par1 = parconvert(wparm1);
      par2 = parconvert(wparm2);
      disp([sprintf('%6.3f %6.3f %6.3f %6.3f %6.4f',par1(1),par1(2),...
          par1(3),par2(3),P),...
          sprintf(' - %3d (%10.5f) %10.5f',iter,errb,dist)]); 
   end
   
   if (size(X,2)>T)
       X = 1:T;
   end
   
   %************* evaluate error in 1st 50ms of fit
   pfit = P*Fisi1(X) + (1-P)*Fisi2(X);
   err = sum((fitisi(X) - pfit) .^ 2);
   tot = sum( (fitisi(X) - mean(fitisi(X)) ).^ 2);
   varexp = ((tot-err)/tot)*100;
   
   pfit = P*Fisi1(X) + (1-P)*Fisi2(X);
   err = sum((tisi(X) - pfit) .^ 2);
   tot = sum( (tisi(X) - mean(tisi(X)) ).^ 2);
   tvarexp = ((tot-err)/tot)*100;
  
   %*****************************************
   pfit = (P*Fisi1) + ((1-P)*Fisi2); 
   ito.worst_loglike = 0;
   ito.best_loglike = 0;
   ito.real_loglike = 0;
   for ii = 1:T
       if (fitisi(ii)>0)
          ito.best_loglike = ito.best_loglike - fitisi(ii) * log(fitisi(ii));
          val = (1/T);
          ito.worst_loglike = ito.worst_loglike - fitisi(ii) * log( val );
          val = pfit(ii);
          ito.real_loglike = ito.real_loglike - fitisi(ii) * log( val );
       end
   end
   
   ito.wparm1 = wparm1;
   ito.wparm2 = wparm2;
   ito.wpar1 = parconvert(wparm1);
   ito.wpar2 = parconvert(wparm2);
   ito.op = OP;
   ito.p = P;
   ito.varexp = varexp;
   ito.tvarexp = tvarexp;
   ito.pfit = pfit;
   ito.errb = errb;
   
   if (plotit >= 1)
       %*********************************************
       subplot(3,3,4);
       plot(XX,fitisi(XX),'k.',XX,(P*Fisi1(XX) + (1-P)*Fisi2(XX)),'k-');
       xlabel('Time (ms)');
       ylabel('Frequency');
       title('Single Fit of ISI');
       subplot(3,3,5);
       plot(X,fitisi(X),'k.',X,(P*Fisi1(X) + (1-P)*Fisi2(X)),'k-');
       xlabel('Time (ms)');
       ylabel('Frequency');
       title(sprintf('Single Fit in 50ms, VarExp: %10.5f',varexp));
       subplot(3,3,6);
       semilogy(X,fitisi(X),'k.',X,(P*Fisi1(X) + (1-P)*Fisi2(X)),'k-');
       xlabel('Time (ms)');
       ylabel('Log Frequency');
       title(sprintf('Single Pars: %6.3f %6.3f %6.3f',par1(1),par1(2),par1(3)));
   end  

 return;


%*************************************************
function ito = dual_fit(fitisi,tisi,awref,bwref,pref,plotit,plotnum,capone,withref)
%******** fit recovery function and poisson isi dist
 global maxit Memb MembR Isi Fisi1 Fisi2 AFisi AFisi1 AFisi2 T WREF TP;

   %******************************************
   crit = 0.0000001;
   % maxit = 300;
   opts = optimset('fminsearch'); %get default values
   optimset(opts,'MaxFunEvals',(10000*6));  % default was 200*6
   
   %******************************************
   disp(sprintf('\n Running gradient descent from best guess'));
    
   %***************************************
   XX = 1:size(fitisi,2);
   Xt = (1/TP) * XX;
   x = find( Xt < 50 );  % find only first 50 ms to plot
   X = 1:size(x,2);
   
   %*******************************************************
   dist = 10000000;
   wrefa = [awref(1),awref(2)];
   wrefb = [awref(1),awref(2)];  % slower start point
   wparm1 = [wrefa,awref(3)];    % single time scale
   wparm2 = [wrefb,bwref(3)];    % don't fit longer time scale
   Fisi1 = isifun(wparm1);
   Fisi2 = isifun(wparm2);
   Isi = fitisi(1:T);
   Isi = Isi / sum(sum( Isi ));
   OP = log( ((1-pref)/pref) );
   P = 1/(1+exp(OP));
   [Memb,MembR] = membifun(wparm1,wparm2,P);  % initial comp of member func
   oerr = isicost3(OP);  % initial estimate of error
   iter = 1;
   while (dist>crit) & (iter<maxit)
       
      %********** fit wparm1 ***************
      WREF = [wrefa,wparm1(3),wparm2(3),P,wrefb];
      wpar1 = wparm1(3);
      if (capone == 0)
        [wpar1,err1] = fminsearch( @isicost1a, wpar1, opts);
      end
      WREF(3) = wpar1;
      wpar2 = wparm2(3);
      WREF(1) = wrefb(1);
      WREF(2) = wrefb(2);
      [wpar2,err3] = fminsearch( @isicost1b, wpar2, opts);
      WREF(4) = wpar2;
      WREF(1) = wrefa(1);
      WREF(2) = wrefa(2);
      
      wparm1(1) = wrefa(1);
      wparm1(2) = wrefa(2);
      wparm1(3) = wpar1;
      wparm2(1) = wrefb(1);
      wparm2(2) = wrefb(2);
      wparm2(3) = wpar2;
      %*************************************
      Fisi1 = isifun(wparm1);
      Fisi2 = isifun(wparm2);
      %********** fit membership ***********
      [OP,errb] = fminsearch( @isicost3, OP, opts);
      errb = isicost3(OP);
      P = 1/(1+exp(OP));
      WREF(5) = P;
      [Memb,MembR] = membifun(wparm1,wparm2,P);
      %*************************************
      
      %******** fit refractory only after others
      if (withref == 1)
        [wrefa,err3] = fminsearch( @isicost1c, wrefa, opts);
        WREF(1) = wrefa(1);
        WREF(2) = wrefa(2);
        wparm1(1) = wrefa(1);
        wparm1(2) = wrefa(2);
        WREF(6) = wrefa(1);
        WREF(7) = wrefa(2);
        wparm2(1) = wrefa(1);
        wparm2(2) = wrefa(2);
        wrefb = wrefa;
        
        [Memb,MembR] = membifun(wparm1,wparm2,P);
        
        Fisi1 = isifun(wparm1);
        Fisi2 = isifun(wparm2);
        errb = isicost3(OP);
      end
      %******* refractory trades off with Tau1, so its a problem
      %******* try to reduce it by fitting Tau1 first *********
      %*********************************************************
      
      dist = (abs(oerr - errb)/abs(oerr)) * 100;
      iter = iter + 1;
      oerr = errb;
      %************************************
      par1 = parconvert(wparm1);
      par2 = parconvert(wparm2);
      disp([sprintf('%6.3f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f',par1(1),par1(2),...
          par1(3),par2(3),P,par2(1),par2(2)),...
          sprintf(' - %3d (%10.5f) %10.5f',iter,errb,dist)]); 
   end
   
   pfit = P*Fisi1(X) + (1-P)*Fisi2(X);
   err = sum((fitisi(X) - pfit) .^ 2);
   tot = sum( (fitisi(X) - mean(fitisi(X)) ).^ 2);
   varexp = ((tot-err)/tot)*100;
   
   pfit = P*Fisi1(X) + (1-P)*Fisi2(X);
   err = sum((tisi(X) - pfit) .^ 2);
   tot = sum( (tisi(X) - mean(tisi(X)) ).^ 2);
   tvarexp = ((tot-err)/tot)*100;
   
   pfit = (P*Fisi1) + ((1-P)*Fisi2);
   ito.worst_loglike = 0;
   ito.best_loglike = 0;
   ito.real_loglike = 0;
   for ii = 1:T
       if (fitisi(ii)>0)
          ito.best_loglike = ito.best_loglike - fitisi(ii) * log(fitisi(ii));
          val = (1/T);
          ito.worst_loglike = ito.worst_loglike - fitisi(ii) * log( val );
          val = pfit(ii);
          ito.real_loglike = ito.real_loglike - fitisi(ii) * log( val );
       end
   end

   ito.wparm1 = wparm1;
   ito.wparm2 = wparm2;
   ito.wpar1 = parconvert(wparm1);
   ito.wpar2 = parconvert(wparm2);
   ito.op = OP;
   ito.p = P;
   ito.varexp = varexp;
   ito.tvarexp = tvarexp;
   ito.pfit = pfit;
   ito.errb = errb;
   
   if (plotit >= 1)
       %*********************************************
       subplot(3,3,7);
       plot(XX,fitisi(XX),'k.',XX,(P*Fisi1(XX) + (1-P)*Fisi2(XX)),'k-');
       xlabel('Time (ms)');
       ylabel('Frequency');
       title('Dual Fit ');
       subplot(3,3,8);
       plot(X,fitisi(X),'k.',X,(P*Fisi1(X) + (1-P)*Fisi2(X)),'k-');
       xlabel('Time (ms)');
       ylabel('Frequency');
       title(sprintf('Dual Fit in 50ms, VarExp: %10.5f',varexp));
       subplot(3,3,9);
       semilogy(X,fitisi(X),'k.',X,(P*Fisi1(X) + (1-P)*Fisi2(X)),'k-');
       xlabel('Time (ms)');
       ylabel('Log Frequency');
       title(sprintf('Dual Pars:%6.3f %6.3f %6.3f %6.3f P=%6.4f',...
                    par1(1),par1(2),par1(3),par2(3),P));
   end  
  
 return;



 %*************** member function computation ***************
 function [membo,membor] = membifun(wparm1,wparm2,P)
   global Isi T TP;
 
   a = isifun(wparm1);
   b = isifun(wparm2);
   
   membo = zeros(1,T);
   for i = 1:T
      membo(1,i) = (P*a(i))/((P*a(i))+((1-P)*b(i)));
   end
   membor = membo;
   
 return
 
 
 %*************** member function computation ***************
 function [membo,membor] = membifunz(wparm1,wparm2,P)
   global Isi T;
   membo = ones(1,T);
   membor = ones(1,T);
   
 return
 

 %******************* limit parameter range ***********************
 function apara = parconvert(para)
     apara = zeros(1,size(para,2));
     %********** constrain recovery to realistic range (Herz, 2006)
     
     for i = 1:size(para,2)
         apara(1,i) = 0.1+exp(para(1,i));   %just bound above zero
     end
        
 return;
 
 %**************** isi parametric function fitting ***********
 function err = isicost1a(par)
   global T Isi Memb WREF;
   
   fisi = isifun([WREF(1),WREF(2),par]);
   
   err = 0;
   minlog = -500;
   for i = 1:T
      if (Memb(i)>0)
           val = fisi(i);
           if (val>0)
             err = err - (Memb(i)*Isi(i) * log(val));
           else
             err = err - (Memb(i)*Isi(i) * minlog);
           end
      end
   end
   
 return;
 
 %**************** isi parametric function fitting ***********
 function err = isicost1b(par)
   global T Isi Memb WREF;
   
   fisi = isifun([WREF(1),WREF(2),par]);
   
   err = 0;
   minlog = -500;
   for i = 1:T
         if (Memb(i)<1)
           val = fisi(i);
           if (val>0)
             err = err - ((1-Memb(i))*Isi(i) * log(val));
           else
             err = err - ((1-Memb(i))*Isi(i) * minlog);
           end
         end
   end
   
 return;
 
 %**************** isi parametric function fitting ***********
 function err = isicost1c(par)
   global T Isi Memb MembR WREF;
   
   fisi = isifun([par,WREF(3)]);
   fisi2 = isifun([par,WREF(4)]);
   P = WREF(5);
   
   err = 0;
   minlog = -500;
   for i = 1:T
         val = P*fisi(i) + (1-P)*fisi2(i);
         if (val>0)
           err = err - (Isi(i) * log(val));
         else
           err = err - (Isi(i) * minlog);
         end
   end
   
 return;
 
 %**************** isi parametric function fitting ***********
 function err = isicost1ca(par)
   global T Isi Memb MembR WREF;
   
   fisi = isifun([par,WREF(3)]);
   fisi2 = isifun([WREF(6),WREF(7),WREF(4)]);
   P = WREF(5);
   
   err = 0;
   minlog = -500;
   for i = 1:T
         val = P*fisi(i) + (1-P)*fisi2(i);
         if (val>0)
           err = err - (Isi(i) * log(val));
         else
           err = err - (Isi(i) * minlog);
         end
   end
   
 return;
 
 %**************** isi parametric function fitting ***********
 function err = isicost1cb(par)
   global T Isi Memb MembR WREF;
   
   fisi = isifun([WREF(1),WREF(2),WREF(3)]);
   fisi2 = isifun([par,WREF(4)]);
   P = WREF(5);
   
   err = 0;
   minlog = -500;
   for i = 1:T
         val = P*fisi(i) + (1-P)*fisi2(i);
         if (val>0)
           err = err - (Isi(i) * log(val));
         else
           err = err - (Isi(i) * minlog);
         end
   end
   
 return;
 
 
 %**************** isi parametric function fitting ***********
 function err = isicost3(op)
   global T Hist Isi Fisi1 Fisi2 AFisi AFisi1 AFisi2 Memb;
   
   AFisi1 = zeros(1,T);
   AFisi2 = zeros(1,T);
   AFisi = zeros(1,T);
   
   P = 1/(1+exp(op));   % constrain to interval 0 to 1
   
   err = 0;
   minlog = -500;
   for i = 1:T
         AFisi1(i) = AFisi1(i) + Isi(i)*Memb(i);
         AFisi2(i) = AFisi2(i) + Isi(i)*(1-Memb(i));
         val = P * Fisi1(i) + (1-P) * Fisi2(i);
        
         AFisi(i) = AFisi(i) + Isi(i)*val;
         
         if (val>0)
           err = err - (Isi(i) * log(val));
         else
           err = err - (Isi(i) * minlog);
         end
   end
   AFisi1 = AFisi1 / sum(sum(AFisi1));
   if (sum(sum(AFisi2))>0)
      AFisi2 = AFisi2 / sum(sum(AFisi2));
   end
   AFisi = AFisi / sum(sum(AFisi));

   err = err + (0.000001 * (op^2));
   
 return;

 
 %********** convolution of two poisson processes
 function ispret = isifun(par)
    global T TP;
    
    MYTP = TP;  % time divisions used to make ISI distribution
    
    UT = 1000*MYTP;  % extra long to insure no normalization cutoff
    
    ispret = zeros(1,T);  % will be converted into ISIT scale at end
    
    ispb = zeros(1,UT);
    isp = zeros(1,UT);
    paro = parconvert(par);
   
    ispa = playerlang(paro(1),paro(2),MYTP);
    MT = size(ispa,2);
    
    %********* there is a fast way to compute this convolution
    %********* because of the exponentials involvement
    papo = (1000/paro(3));
    if (papo==0)
          isp = ispa;  % exponent becomes delta in limit
    else
          isp(1) = ispa(1)*1;  % first exponent is 1 always, in limit, defaults to ispa
          pb = exp( -((1/MYTP)/papo));
          for i = 2:UT   % freakin fast compared to regular way
             if (i<=MT) 
               isp(i) = isp(i-1)*pb + ispa(i)*pb;
             else
               isp(i) = isp(i-1)*pb;
             end
             if (i>MT) & (isp(i) == 0)
                 break;  % dont continue if down to zero already
             end
          end
    end     
    
    %*****************************************************************
    isp = isp / sum(sum( isp )); 
    
    ispret = isp;  % should match identical timing
    
 return;
 
 
%***************************************************************
function gamma = playerlang(tauref,erpow,MYTP)
  global T TP;

  UT = 2000;
  
  tr = tauref;   % mean refractory is tau1 thus, but erlang distribution
  %********** I'm using the Michelson-Menten Recovery Function from Herz
  %********** which has a more constrained form than the erlang *******
  igamma = zeros(1,UT);  %let size of gamma be variable (as kernel size could be very small)
  for i = 1:UT
        ta = (i/MYTP);
        val = ((-ta+tauref)/erpow);  
        w = 1/(1+exp(val));
        dw = w*(1-w);
        igamma(i) = dw; 
        if (ta>tauref) & (dw==0)
          break;
        end
  end
  gamma = igamma(1:i);
  
  zsum = sum(sum(gamma));
  if (zsum>0)
     gamma = gamma / zsum;
  else
     gamma = 1;  %default to delta function in limit
  end

return;


 %*************** compute raw isi distributions ***************
 function [isibin,isicnt,redflag] = comp_isifreq(spikearray)
   global TP;  % time divisions should be 0.5 ms
   pad = 0;  
   
   N = size(spikearray,2);
   
   allspikes = spikearray{1};
   M = floor( (size(allspikes,2)-pad)*3/4);
   pstart = pad+1;
   pend = size(allspikes,2);
  
   
   redflag = 0;
   %****************************** for norm scale, use 0.5ms divisions
   %normisiscale = (1:floor((2*M)/TP)) * 0.5;
   %M2 = size(normisiscale,2);
   M2 = M;
   isibin = zeros(N,M2);
   isicnt = zeros(N,M2);  
   
   %**************** find the ISI distributions ********************
   %********* (also, correct for there being less samples of longer
   %*********  ISI durations due to finite size of data interval by
   %*********  dividing by the total number of possible samples at
   %*********  each ISI, for example, for 800ms, you can sample a
   %*********  1ms ISI 800 times, but a 400ms isi only 400 times).
   for zkk = 1:N
     zisicnt = zeros(1,M2);
     zisibin = zeros(1,M2);
     allspikes = spikearray{zkk};
     %**********************************************
     for ii = 1:size(allspikes,1)  %ii is the trial number, 1-160 for all spikes
       if (mod(ii,10) == 0)
         disp(sprintf('Computing ISI distribution, trial: %d',ii));
       end
       jj = pstart;
       while (jj<pend)
           if (allspikes(ii,jj)==1)
               hito = jj;
               for kk = (jj+1):pend
                 if ((kk-jj)<=M)
                    it = ((kk-jj)/TP); 
                    ikk = (kk-jj);
                    %*********************************************  
                    zisicnt(1,ikk) = zisicnt(1,ikk) + 1;
                 end
               end
               for kk = (jj+1):pend
                  if (allspikes(ii,kk)==1)
                      hito = kk;
                      break;
                  end
               end
               if (hito > jj)  %found next spike
                  if ((hito-jj)<=M)
                    it = ((hito-jj)/TP); 
                    ikk = (hito-jj);
                    zisibin(1,ikk) = zisibin(1,ikk) + 1;
                  end
                  jj = hito;
               else
                  jj = jj + 1;
               end
           else
               jj = jj + 1;
           end
        end %jj updated
     end
     
     %***************** normalize by test frequency ***********
     for ii = 1:size(zisibin,2)
       if (zisicnt(1,ii)>0)
         zisibin(1,ii) = zisibin(1,ii) / zisicnt(1,ii);
       else
         sprintf('ii has zero count %d');
       end
     end
     if (sum(sum(zisibin(1,:))) == 0)
        zisibin(1,:) = 0;
        redflag = zkk;
     else
        zisibin(1,:) = zisibin(1,:) / sum(sum(zisibin(1,:)));
     end
     %*******************************************************
     isibin(zkk,:) = zisibin; 
     isicnt(zkk,:) = zisicnt;   
  end  %zkk looped
   
 return;
 
 