function [SCORE] = change_detection(X, n, k, alpha)

SCORE=[];

WIN = sliding_window(X, k, 1);
nSamples = size(WIN,2);
t = n +1;

h = waitbar(0,'Please wait...');
while(t + n -1 <= nSamples)
    waitbar((t + n -1) / nSamples);
    
    YTest = WIN(:, t : n + t -1 );
    YRef = WIN(:, t-n: t-1);
    [pe]=RelULSIF(YRef,YTest,[],[],alpha);
    
    SCORE = [SCORE, pe];
    t = t + 1;
end
close(h)
end

