function RESP = SPK_resp_profile(spks, SRT, Rew, do_clip, do_plot)
% quatify the response profile of a neuron.
%
% DESCRIPTION
% Get some numbers that quatify the response profile of a neuron> The input is
% assumed to be a trial matrix with spike times aligned to stimulus onset.
%
% SYNTAX
%   RESP = SPK_resp_profile(spk, movev, do_plot)
%
%   Input:
%         <spk>        spike time or spike object
%
%         <SRT>        second event time (e.g. saccade) to align the responses
%
%         <do_clip>    clip spikes after saccade from stimulus aligned trials
%
%         <do_plot>    plot the results
%
%
% REFERENCES
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% wolf zinke, 20.2.2014

% ____________________________________________________________________________ %
%% define time windows for the different response periods
Spont   = [-250;   0];

VisTran = [  50; 100];
VisSust = [ 100; 150];

MovSust = [-125; -75];
MovTran = [ -50;   0];
MovSacc = [   0;  50];
MovPost = [  50; 100];

PreRew  = [ -50;   0];
PostRew = [  50; 100];

uthr = 90; % percentile threshold based on spont activity to be accountes as responses modulation
lthr = 10;

%% check data input
if(~exist('Rew','var'))
    Rew = [];
end

if(~exist('do_clip','var') || isempty(do_clip))
    do_clip = 0;
end

if(~exist('do_plot','var') || isempty(do_plot))
    do_plot = 0;
end

Ntrials = size(spks,1);

if(size(spks,2) < 4)
    warning('Not enough trials to do something useful here!');
    RESP = [];
    return;
end

if(Ntrials ~= length(SRT))
    error('Number of SRT values has to be the same as tiral number!');
end

sacc_spks = bsxfun(@minus, spks, SRT(:));

% ____________________________________________________________________________ %
%% get responses for the different epochs
[RESP.Spont_mean, RESP.Spont_var, RESP.Spont_ste, RESP.Spont_fano, RESP.Spont_rate] = get_meanresp(spks, Spont);

% clip spikes after saccade (subtract 10 to account for minor synaptic latency)
if(do_clip == 1)
    [RESP.VisTran_mean, RESP.VisTran_var, RESP.VisTran_ste, RESP.VisTran_fano, RESP.VisTran_rate] = get_meanresp(spks, VisTran, SRT(:)-10);
    [RESP.VisSust_mean, RESP.VisSust_var, RESP.VisSust_ste, RESP.VisSust_fano, RESP.VisSust_rate] = get_meanresp(spks, VisSust, SRT(:)-10);
else
    [RESP.VisTran_mean, RESP.VisTran_var, RESP.VisTran_ste, RESP.VisTran_fano, RESP.VisTran_rate] = get_meanresp(spks, VisTran);
    [RESP.VisSust_mean, RESP.VisSust_var, RESP.VisSust_ste, RESP.VisSust_fano, RESP.VisSust_rate] = get_meanresp(spks, VisSust);
end

[RESP.MovSust_mean, RESP.MovSust_var, RESP.MovSust_ste, RESP.MovSust_fano, RESP.MovSust_rate] = get_meanresp(sacc_spks, MovSust);
[RESP.MovTran_mean, RESP.MovTran_var, RESP.MovTran_ste, RESP.MovTran_fano, RESP.MovTran_rate] = get_meanresp(sacc_spks, MovTran);
[RESP.MovSacc_mean, RESP.MovSacc_var, RESP.MovSacc_ste, RESP.MovSacc_fano, RESP.MovSacc_rate] = get_meanresp(sacc_spks, MovSacc);
[RESP.MovPost_mean, RESP.MovPost_var, RESP.MovPost_ste, RESP.MovPost_fano, RESP.MovPost_rate] = get_meanresp(sacc_spks, MovPost);

if(~isempty(Rew))
    [RESP.PreRew_mean,  RESP.PreRew_var,  RESP.PreRew_ste,  RESP.PreRew_fano,  RESP.PreRew_rate]  = get_meanresp(spks,  PreRew);
    [RESP.PostRew_mean, RESP.PostRew_var, RESP.PostRew_ste, RESP.PostRew_fano, RESP.PostRew_rate] = get_meanresp(spks, PostRew);    
end

RESP.VisTran_norm = RESP.VisTran_mean / RESP.Spont_mean;
RESP.VisSust_norm = RESP.VisSust_mean / RESP.Spont_mean;
RESP.MovSust_norm = RESP.MovSust_mean / RESP.Spont_mean;
RESP.MovTran_norm = RESP.MovTran_mean / RESP.Spont_mean;
RESP.MovSacc_norm = RESP.MovSacc_mean / RESP.Spont_mean;
RESP.MovPost_norm = RESP.MovPost_mean / RESP.Spont_mean;

if(~isempty(Rew))
    RESP.PreRew_norm  = RESP.PreRew_mean  / RESP.Spont_mean;
    RESP.PostRew_norm = RESP.PostRew_mean / RESP.Spont_mean;
end

% get response modulation and identify differences to spont activity
% resp_thr = RESP.Spont_mean + 2 * sqrt(RESP.Spont_var);
resp_thr = RESP.Spont_mean + 2 * RESP.Spont_ste;

if(RESP.VisTran_mean > resp_thr || RESP.VisSust_mean > resp_thr)
    RESP.VisResp = 1;
else
    RESP.VisResp = 0;
end

if(RESP.MovTran_mean > resp_thr)
    RESP.MovResp = 1;
else
    RESP.MovResp = 0;
end

if(~isempty(Rew))
    if(RESP.PostRew_mean > RESP.PreRew_mean + 2 * RESP.PreRew_ste)
        RESP.RewResp = 1;
    else
        RESP.RewResp = 0;
    end
end

if(RESP.VisResp == 1 && RESP.MovResp == 1)
    RESP.resptype = 'vismov';
elseif(RESP.VisResp == 1)
    RESP.resptype = 'vis';
elseif(RESP.MovResp == 1)
    RESP.resptype = 'mov';
else
    RESP.resptype = 'none';
end

% do the same based percentiles to get also a suppression
upperthr = prctile(RESP.Spont_rate, uthr);
lowerthr = prctile(RESP.Spont_rate, lthr);

if(median(RESP.VisTran_rate) > upperthr )
    RESP.VisTran_Mod = 1;
elseif(median(RESP.VisTran_rate) < lowerthr)
    RESP.VisTran_Mod = -1;
else
    RESP.VisTran_Mod = -0;
end

if(median(RESP.VisSust_rate) > upperthr )
    RESP.VisSust_Mod = 1;
elseif(median(RESP.VisSust_rate) < lowerthr)
    RESP.VisSust_Mod = -1;
else
    RESP.VisSust_Mod = -0;
end

if(median(RESP.MovTran_rate) > upperthr )
    RESP.MovTran_Mod = 1;
elseif(median(RESP.MovTran_rate) < lowerthr)
    RESP.MovTran_Mod = -1;
else
    RESP.MovTran_Mod = -0;
end

if(median(RESP.MovSacc_rate) > upperthr )
    RESP.MovSacc_Mod = 1;
elseif(median(RESP.MovSacc_rate) < lowerthr)
    RESP.MovSacc_Mod = -1;
else
    RESP.MovSacc_Mod = -0;
end

if(median(RESP.MovPost_rate) > upperthr )
    RESP.MovPost_Mod = 1;
elseif(median(RESP.MovPost_rate) < lowerthr)
    RESP.MovPost_Mod = -1;
else
    RESP.MovPost_Mod = -0;
end

% ____________________________________________________________________________ %
%% calculate ratios for the responses

% visual transiency
RESP.VisTrans_IDX   = get_contrast(RESP.VisTran_mean, RESP.VisSust_mean);

% saccade transiency
RESP.MovTrans_IDX   = get_contrast(RESP.MovTran_mean, RESP.MovSust_mean);

% saccade periods
RESP.MovPrePost_IDX  = get_contrast(RESP.MovTran_mean, RESP.MovPost_mean);
RESP.MovPreSacc_IDX  = get_contrast(RESP.MovTran_mean, RESP.MovSacc_mean);
RESP.MovSaccPost_IDX = get_contrast(RESP.MovSacc_mean, RESP.MovPost_mean);

% mov/vis contrast
RESP.VisMov_IDX   = get_contrast(RESP.VisTran_mean, RESP.MovTran_mean);
RESP.VisSpont_IDX = get_contrast(RESP.VisTran_mean, RESP.Spont_mean);
RESP.MovSpont_IDX = get_contrast(RESP.MovTran_mean, RESP.Spont_mean);

% ____________________________________________________________________________ %
%% get p-values for the comparisons as well
% use a paired Wilcoxon signed rank test, because samples are not independent
% from each other but from the same trial.

% visual transiency
RESP.VisTrans_P    = signrank(RESP.VisTran_rate, RESP.VisSust_rate);

% saccade transiency
RESP.MovTrans_P    = signrank(RESP.MovTran_rate, RESP.MovSust_rate);

% saccade periods
RESP.MovPrePost_P  = signrank(RESP.MovTran_rate, RESP.MovPost_rate);
RESP.MovPreSacc_P  = signrank(RESP.MovTran_rate, RESP.MovSacc_rate);
RESP.MovSaccPost_P = signrank(RESP.MovSacc_rate, RESP.MovPost_rate);

% mov/vis contrast
RESP.VisMov_P      = signrank(RESP.VisTran_rate, RESP.MovTran_rate);

% response modulation
RESP.VisSpont_P    = signrank(RESP.VisTran_rate, RESP.Spont_rate);
RESP.MovSpont_P    = signrank(RESP.MovTran_rate, RESP.Spont_rate);

% ____________________________________________________________________________ %
%% plot a summary representation
if(do_plot == 1)
    figure;

% ____________________________________________________________________________ %
    subplot(3,4,1);
    wz_spk_plot_hist(spks, [-50, 150], [],0,10,0,100,1/3, SRT(:)-10);
    title('Stimulus', 'FontSize', 12, 'Interpreter', 'none');

% ____________________________________________________________________________ %
    subplot(3,4,2);
    wz_spk_plot_hist(sacc_spks, [-150, 50], [],0,10,0,100,1/3);
    title('Saccade', 'FontSize', 12, 'Interpreter', 'none');

% ____________________________________________________________________________ %
    subplot(3,4,3);

    mnrt = [RESP.Spont_mean;   RESP.VisTran_mean; RESP.VisSust_mean; RESP.MovTran_mean; ...
            RESP.MovSust_mean; RESP.MovSacc_mean; RESP.MovPost_mean];

    bar(mnrt);
    set(gca,'XTickLabel',{'Spont'; 'VisTran'; 'VisSust'; 'MovTran'; 'MovSust'; 'MovSacc'; 'MovPost'})
    rotateXLabels( gca(), 45);
    title('Mean Rate [spk/]', 'FontSize', 12, 'Interpreter', 'none');


% ____________________________________________________________________________ %
    subplot(3,4,4);

    fano = [RESP.Spont_fano;   RESP.VisTran_fano; RESP.VisSust_fano; RESP.MovTran_fano; ...
            RESP.MovSust_fano; RESP.MovSacc_fano; RESP.MovPost_fano];

    bar(fano);
    set(gca,'XTickLabel',{'Spont'; 'VisTran'; 'VisSust'; 'MovTran'; 'MovSust'; 'MovSacc'; 'MovPost'})
    rotateXLabels( gca(), 45);
    title('Fano Factor', 'FontSize', 12, 'Interpreter', 'none');

% ____________________________________________________________________________ %
    subplot(3,4,5);
    plot_scatter(RESP.VisTran_rate, RESP.VisSust_rate, 'Vis Transiency', 'VisTran','VisSust');

% ____________________________________________________________________________ %
    subplot(3,4,9);
    plot_scatter(RESP.MovTran_rate, RESP.MovSust_rate, 'Mov Transiency', 'MovTran','MovSust');

% ____________________________________________________________________________ %
    subplot(3,4,6);
    plot_scatter(RESP.VisTran_rate, RESP.MovTran_rate, 'Vis/Mov', 'VisTran','MovTran');

% ____________________________________________________________________________ %
    subplot(3,4,7);
    plot_scatter(RESP.Spont_rate, RESP.VisTran_rate, 'Vis', 'Spont','VisTran');

% ____________________________________________________________________________ %
    subplot(3,4,8);
    plot_scatter(RESP.Spont_rate, RESP.MovTran_rate, 'Mov', 'Spont','MovTran');

% ____________________________________________________________________________ %
    subplot(3,4,10);
    plot_scatter(RESP.MovTran_rate, RESP.MovPost_rate, 'Pre/Post', 'Pre Saccade','Post Saccade');

% ____________________________________________________________________________ %
    subplot(3,4,11);
    plot_scatter(RESP.MovTran_rate, RESP.MovSacc_rate, 'Pre/Sacc', 'Pre Saccade','Saccade');

% ____________________________________________________________________________ %
    subplot(3,4,12);
    plot_scatter(RESP.MovSacc_rate, RESP.MovPost_rate, 'Sacc/Post', 'Saccade','Post Saccade');

end

% ____________________________________________________________________________ %
%% helper functions
function ctr = get_contrast(a,b)
    ctr = (a - b) / (a + b);


function [meanrt, varrt, stert, fano, trialrate] = get_meanresp(spktm, timwin, clip)

    if(exist('clip','var') && ~isempty(clip))
        pos = spktm >= timwin(1) & bsxfun(@le, spktm, clip(:));
        spkcnt    = sum(pos,2);
        trialrate = 1000 .* spkcnt ./ (clip(:) - timwin(1) +1);
    else
        pos = spktm >= timwin(1)  & spktm <= timwin(2);
        spkcnt    = sum(pos,2);
        trialrate = 1000 .* spkcnt ./ (diff(timwin)+1);
    end

    meanrt = mean(trialrate);
    varrt  =  var(trialrate);
    stdrt  =  std(trialrate);
    stert  =  stdrt / sqrt(length(trialrate));
    fano   = varrt/meanrt;


function plot_scatter(x,y,ttl, xlab, ylab)

    minVal = min([x(:);y(:)]);
    maxVal = max([x(:);y(:)]);

    medX   = prctile(x,[50,24,75]);
    medY   = prctile(y,[50,24,75]);

    st = 0.05 * (maxVal - minVal);

    plot([minVal maxVal],[minVal maxVal],'--r','LineWidth',1);
    hold on;

    plot(x, y, 'ok','MarkerFaceColor','k', 'MarkerSize',4, 'MarkerEdgeColor','k');

    plot(medX(1), medY(1), 'og','MarkerFaceColor','g', 'MarkerSize',8, 'MarkerEdgeColor','g');
    plot(medX(2:3),medY([1,1]),'-g','LineWidth', 2);
    plot(medX([1,1]),medY(2:3),'-g','LineWidth', 2);

    xlim([minVal-st, maxVal+st]);
    ylim([minVal-st, maxVal+st]);

    axis square

    title(ttl, 'FontSize', 12, 'Interpreter', 'none');
    xlabel(xlab, 'FontSize', 10, 'Interpreter', 'none');
    ylabel(ylab, 'FontSize', 10, 'Interpreter', 'none');


