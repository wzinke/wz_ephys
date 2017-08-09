clear all;
addpath RULSIF

figure
load logwell.mat

% choice of alpha
% better estimation,  but look messier (since divergence is upperbounded)
% alpha = .1;

% divergence is not bounded, less accurate, but good for visualization
alpha = .0;

n = 50;
k = 10;

score1 = change_detection(y,n,k,alpha);
score2 = change_detection(y(:,end:-1:1),n,k,alpha);

subplot(2,1,1);
plot(y, 'b-', 'linewidth',2);
axis([-inf,size(y,2),-inf,inf])
title('Original Signal')

subplot(2,1,2);
score2 = score2(end:-1:1);

% 2*n+k-2 is the size of the "buffer zone".
plot([zeros(1,2*n-2+k),score1 + score2], 'r-', 'linewidth',2);
axis([-inf,size(y,2),-inf,inf])
title('Change-Point Score')
