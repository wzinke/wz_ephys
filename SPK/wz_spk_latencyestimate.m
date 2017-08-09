function wz_spk_latencyestimate(spk, timwin, nBoot, limits)
% estimate response latency based on a description given by Gawne et al. (1996)
%    - convolution of individual spikes with a gaussian kernel
%    - averaging rate estimates over all trials
%    - determination of half-peak time
%    - if half-peak is at least twice of spontaneous activity latency is valid
%    -> added to original description:
%        - bootstrap validation of the latency
%
% wolf zinke, 26.10.2004

est_level  = 'poiss';   % What level is chosen to asses whether response is present
level_fact = 0.05;      % how many times does it need to be higher than the level

%% check data input
if(~isstruct(spk))
    if(length(unique(spk(~isnan(spk)))) > 2)
        spk = wz_get_SPKobj(spk);
    end
end

if(~exist('timwin','var') || isempty(timwin))
    timwin = [0, spk.timewindow(2)];
end

if(~exist('nBoot','var') ||isempty(nBoot) == 1)
    nBoot = 1;
end

if(~exist('thresh_val','var') || isempty(thresh_val) == 0)
    thresh_val = [];
end

% Time window within a valid latency is expected
if(~exist('limits','var') || isempty(limits) == 1)
    limits = [20;150];
end

spk = wz_spk_density(spk, 'gauss', 5);

    [resp_estimate, ste_trial, var_vals, sgl_vals] = WZ_trial_rates(spiketrain, krnl_width(j),krnl_type,use_half);

    % define characteristics of spontaneous activity
    spont_bin = find(bintimes < limits(1));
    mean_spont = mean(resp_estimate(spont_bin));
    spont_var_level = mean(var_vals(spont_bin));
    spont_std_level = mean((var_vals(spont_bin)).^0.5);

if(isempty(thresh_val) == 1)
    switch est_level
        case 'spont'
            use_level = level_fact * mean_spont;

        case 'var'
            use_level = level_fact * spont_var_level;

        case 'std'
            use_level = level_fact * spont_std_level;

        case 'mad'
            use_level = level_fact * spont_std_level;
            
        case 'poiss'
            if(level_fact >=1)
                level_fact = 0.05;
            end
            [~, ci_vals] = poissfit(resp_estimate(spont_bin),level_fact);
            use_level = ci_vals(2);

        otherwise
            error(['Chosen level <',est_level,'> not known!']);
    end
else
    use_level = thresh_val;
end



figure()
wz_spk_plot_density(spk,[-200,400])
