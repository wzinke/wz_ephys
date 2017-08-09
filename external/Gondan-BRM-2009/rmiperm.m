function [p, tmax, tcrit] = rmiperm(fname, onlyobs, q, nperm)
% Permutation test of the race model inequality
%
% fname: filename to be read
%
%     obs cond    rt
%       1    A   234
%       1   AV   263
%       1    V   298
%       1    A   Inf
%
%     obs:  Observer number
%     cond: Condition (A, V, AV)
%     rt:   response time
%     Inf:  Omitted response
%
% onlyobs (optional): restrict analysis to these observers in the analysis
%     default []: use all observers
%
% q: quantiles of overall response time distribution at which RMI should be 
%    tested, default 0.05:0.05:0.30
%
% nperm: number of permutations used to determine tmax distribution
%     default 10001
%
% Tested with Matlab version 7.0.1
%

% defaults if parameters are omitted

if nargin < 1
    fname = 'c:\rtdata\go2004.dat' 
end

if nargin < 2
    onlyobs = [] 
end

if nargin < 3
    q = 0.05:0.05:0.30 
end

if nargin < 4
    nperm = 10001
end

% default return values

p = -1 ;
tcrit = -1 ;
tmax = 'error' ;

% read file

[fid, message] = fopen(fname, 'r') ;
if fid == -1
    disp(message)
    return
end

header = textscan(fid, '%s %s %s', 1) ;
if header{1}{1} ~= 'obs'
    disp('Error reading header.')
    return
end

if header{2}{1} ~= 'cond'
    disp('Error reading header.')
    return
end

if header{3}{1} ~= 'rt'
    disp('Error reading header.')
    return
end

data = textscan(fid, '%f %s %f', 'headerlines', 1) ;

fclose(fid) ;

% column 1: observer no
% column 2: condition (A, V, AV)
% column 3: response time

obs  = data{1} ;
cond = data{2} ;
rt   = data{3} ;

% Extract response times

if isempty(onlyobs)
    onlyobs = min(obs):max(obs) ;
end

d = [] ;

for i = onlyobs
    tA  = rt(obs == i & strcmp('A', cond)) ;
    tV  = rt(obs == i & strcmp('V', cond)) ;
    tAV = rt(obs == i & strcmp('AV', cond)) ;
    tAll = sort(rt(obs == i)) ;

    % Quantiles of overall RTs
    edges = tAll(round(q*length(tAll))) ;

    % Empirical CDFs
    maxrt = 1000 ;                  % adjust if necessary
    fA   = histc(tA, 1:maxrt) ;
    fV   = histc(tV, 1:maxrt) ;
    fAV  = histc(tAV, 1:maxrt) ;

    cdfA  = cumsum(fA)/length(tA) ;
    cdfV  = cumsum(fV)/length(tV) ;
    cdfAV = cumsum(fAV)/length(tAV) ;

    % Race Model Difference
    d = [d ; cdfAV(edges)' - cdfA(edges)' - cdfV(edges)'] ;
end

% Observed tmax statistic

tmax = max(mean(d) ./ std(d)) * sqrt(length(onlyobs)) ;

% Distribution of tmax

tmaxi = zeros(1, nperm) ;
for i = 1:nperm
    dp       = diag(sign(rand(1, length(onlyobs)) - 0.5)) * d ;
    tmaxi(i) = max(mean(dp) ./ std(dp))  * sqrt(length(onlyobs)) ;
end

% 95th quantile of simulated tmax

tmaxi = sort(tmaxi) ;
tcrit = tmaxi(round(0.95*nperm)) ;

% P value

p = sum(tmax < tmaxi)/nperm ;

if nargout < 3
    disp(sprintf('tmax=%.3f tcrit=%.3f P=%.3f', tmax, tcrit, p))
end

return
