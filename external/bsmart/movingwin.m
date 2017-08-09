%function [R1,R2,diffR,ratio]=movingwin(gelb, stim, shif, arcoeff, arnoise);
function [R1,R2,diffR,ratio]=movingwin(gelb,arcoeff, arnoise,WIN1,NPTS1,NTRLS1,NCHN1);
% Input:
%   gelb: T x chans x trials; T - length of the data
%   stim:  absolute time of stimulus, Trials x 1
%   shif:  absolute time of shift, Trials x 1
%   arcoeff: 101 x 1350
%   arnoise: 101 x 225
% Output:
%   R, Ru, Rd: middle, up, down distance
%


global WIN NPTS NTRLS NCHN 
WIN = WIN1;
NPTS = NPTS1;
NTRLS =NTRLS1;
NCHN = NCHN1;
stimmin=0;shfmin=0;shfmax=0;

R1=[];
R2=[];
diffR=[];
ratio=[];
%t = -(stimmin - WIN/2)*5 - 5;   % - (20-8)*5 +5 = -60 - 5 msec 
for rec=0:(NPTS-WIN-shfmax)  % without +1 in matlab
   %t = t + 5;
   idx = rec-abs(shfmin)+1;
   %fprintf('index, t = %d  %d msec\n', idx, t);
   dat1 = zeros(WIN, NCHN, NTRLS);  % monkey data
   for i=1:NTRLS
	% pp = rec+1;  % plus 1 since index in matlab is 1 
	 dat1(:,:,i) = gelb(rec+1:rec+WIN,:,i);
   end

   mar = MAR_make(arcoeff(idx,:), arnoise(idx,:));
   dat2 = ensembding(mar, WIN, NCHN, NTRLS); % generate simulation data
   [R11, R21, diffR1, ratio1] = constchk(dat1, dat2);
   
   R1=[R1;R11];
   R2=[R2;R21];
   diffR = [diffR;diffR1];
   ratio =[ratio;ratio1];
   
end
 


