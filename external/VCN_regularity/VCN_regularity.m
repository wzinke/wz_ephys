function [mu,sigma,CV,rate] = VCN_regularity(t,spikes,type,N_min)
%VCN_REGULARITY Regularity analysis for Ventral Cochlear Nucleus neurons.
%   MU = VCN_REGULARITY(T,SPIKES) returns the mean Inter-Spike Interval
%   (ISI) at the times contained in the vector T for the spike trains 
%   contained in the columns of the array SPIKES. Each column of SPIKES 
%   should consist of an ordered set of spike times (ordering is not 
%   checked) padded with NaNs. MU will contain NaNs for times at which it 
%   is not defined.
%
%   [MU,SIGMA,CV,RATE] = VCN_REGULARITY(...) also returns the standard 
%   deviation SIGMA, the coefficient of variation CV and the mean
%   inter-spike rate RATE.
%
%   VCN_REGULARITY(T,SPIKES,'hybrid') uses a hybrid algorithm to mimic the
%   results of the regularity analysis used by Young et al (1988). ISIs
%   only contribute to their statistics if they are within a certain
%   distance of their commencement, and that distance is controlled by the
%   internal parameter pc_hybrid. VCN_REGULARITY(T,SPIKES,'exact') gives the 
%   default exact calculation.
%
%   VCN_REGULARITY(T,SPIKES,TYPE,N) replaces the outputs with NaNs at times
%   for which fewer than N intervals are in progress.  If N is not
%   specified it is set to one fifth of the number of columns in SPIKES.
%
%References
%   Young, E. D., Robert, J.-M. & Shofner, W. P. "Regularity and
%   latency of units in ventral cochlear nucleus: implications for unit
%   classification and generation of response properties", Journal of
%   Neurophysiology 60, 1-29 (1988).
%
%   Wright, M. C. M., Bleeck, S. & Winter, I. M. "A note on methods of
%   regularity analysis for cochlear nucleus neurons", in preparation for
%   Journal of the Acoustical Society of America (2011).

% Matthew Wright, Institute of Sound & Vibration Research, University of
% Southampton, mcmw@isvr.soton.ac.uk, 31/2/2011.


if nargin < 3
    type = 'exact';
end

N_sweeps = size(spikes,2);

if nargin < 4
    N_min = N_sweeps/5;
end

R = nan(length(t),N_sweeps);

pc_hybrid = 5; % Percentile of all interval lengths to use as 
               % proximity criterion for hybrid method

if strcmpi(type,'hybrid')
    d = sort(diff(spikes(:)));
    d = d(~isnan(d));
    bw = d(ceil(length(d)*pc_hybrid/100));
end

for j = 1:N_sweeps
    for i = 1:length(find(~isnan(spikes(:,j)))) - 1
        if strcmpi(type,'exact')
            R(t >= spikes(i,j) & t < spikes(i+1,j),j) = spikes(i+1,j) - spikes(i,j);
        elseif strcmpi(type,'hybrid')
            R(t >= spikes(i,j) - bw/2 & t < spikes(i,j) + bw/2,j) = spikes(i+1,j) - spikes(i,j);
        else
            error('Unknown type of regularity analysis')
        end
    end
end

mu = nan(size(t));
sigma = mu;
rate = mu;

for i = 1:length(t)
    s = R(i,~isnan(R(i,:)));
    if length(s) > N_min
        mu(i) = mean(s);
        rate(i) = mean(1./s);
        sigma(i) = std(s,1);
    end
end
CV = sigma./mu;

