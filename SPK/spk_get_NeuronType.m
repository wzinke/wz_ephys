function RESP = spk_get_NeuronType(spks, stim, sacc, rew, loc, krnl, kw, max_rate, max_trial, do_plot)
% spk_get_NeuronType - simple classification of a neuron based on response
%       characterisitcs as visual, movement, vismov, and reward related.
%
% DESCRIPTION
%
% SYNTAX
%   TST = wz_spk_getTSdiff()
%
%   Input:  spks
%           stim
%           sacc
%           rew
%           loc
%           krnl
%           kw
%           max_rate
%           max_trial
%
%   Output:
%
%
% REFERENCES
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 24-Sep-2015 by wolf zinke
% $Modified:

% ____________________________________________________________________________ %
%% define time windows for the different response periods (keep response windows at equal size)

PreFix  = [-250,    0];
Spont   = [-250,    0];

VisTran = [  50,  100];
VisSust = [ 100,  150];
VisLate = [ 150,  200];

MovEarl = [-150, -100];
MovSust = [-100,  -50];
MovTran = [ -50,    0];
MovSacc = [   0,   50];
MovPost = [  50,  100];

EarlRew = [-150, -100];
PreRew  = [ -50,    0];
PostRew = [  50,  100];

% define default colors
% main_col = [     0,      0, 1.0000; ...
%             1.0000,      0,      0; ...
%                  0, 1.0000,      0; ...
%             1.0000, 0.1034, 0.7241; ...
%             1.0000, 0.8276,      0];
  

main_col = [228, 26, 28; ...
             55,126,184; ...
             77,175, 74; ...
            152, 78,163; ...
            255,127,  0; ...
            255,255, 51; ...
            166, 86, 40; ...  
            247,129,191]	./ 255;

ptch_col = 1 - ((1-main_col) .* 0.25);

do_alpha = 0;  % use transparency

%% define default arguments
if(~exist('loc','var') || isempty(loc))
    loc = ones(size(spks,1));
end

if(~exist('krnl','var') || isempty(krnl))
    krnl = 'gauss';
end

if(~exist('kw','var') || isempty(kw))
    kw = 10;
end

if(~exist('max_rate','var') || isempty(max_rate))
    get_rate = 1;
    max_rate = 0;
end

if(~exist('max_trial','var') || isempty(max_trial))
    get_trial = 1;
    max_trial = 0;
end

if(~exist('do_plot','var') || isempty(do_plot))
    do_plot = 1;
end
loclst = sort(unique(loc));

sacc_corr = sacc-stim;
sacc_prct = ceil(prctile(sacc_corr,95)/50)*50;

% ____________________________________________________________________________ %
%% get spike density function
 
VisResp = [];
MovResp = [];
RewResp = [];

RESP.VisResp = VisResp; % initialize here, re-do it later with values
RESP.MovResp = MovResp;
RESP.RewResp = RewResp;

RESP.Cond = ensure_row(loclst);

if(any(spks(:) < 0))
    [RESP.PreFix_mean, RESP.PreFix_var, RESP.PreFix_ste, ~, PreFixrates] = get_meanresp(spks, PreFix);
    prct = prctile(PreFixrates,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.PreFix_P01  = prct(1);
    RESP.PreFix_P02  = prct(2);
    RESP.PreFix_P05  = prct(3);
    RESP.PreFix_P10  = prct(4);
    RESP.PreFix_P25  = prct(5);
    RESP.PreFix_median  = prct(6);
    RESP.PreFix_P75  = prct(7);
    RESP.PreFix_P90  = prct(8);
    RESP.PreFix_P95  = prct(9);
    RESP.PreFix_P98  = prct(10);
    RESP.PreFix_P99  = prct(11);
    RESP.PreFix_MAD  = mad(PreFixrates, 1);
    RESP.PreFixrates = PreFixrates;   
end

AllVis = wz_get_SPKobj(spks, [-500,  1000], stim);
[RESP.Base_mean, RESP.Base_var, RESP.Base_ste, ~, Baserates] = get_meanresp(AllVis.spiketimes, Spont);
prct = prctile(Baserates,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
RESP.Base_P01  = prct(1);
RESP.Base_P02  = prct(2);
RESP.Base_P05  = prct(3);
RESP.Base_P10  = prct(4);
RESP.Base_P25  = prct(5);
RESP.Base_median  = prct(6);
RESP.Base_P75  = prct(7);
RESP.Base_P90  = prct(8);
RESP.Base_P95  = prct(9);
RESP.Base_P98  = prct(10);
RESP.Base_P99  = prct(11);
RESP.Base_MAD  = mad(Baserates, 1);
RESP.Base_rate = Baserates;

for(l=1:length(loclst))
    p = loc == loclst(l);
    VisResp = [VisResp, wz_spk_density(wz_get_SPKobj(spks(p,:), [-500,  1000], stim(p)), krnl, kw)];
    MovResp = [MovResp, wz_spk_density(wz_get_SPKobj(spks(p,:), [-1000, 1000], sacc(p)), krnl, kw)];
    RewResp = [RewResp, wz_spk_density(wz_get_SPKobj(spks(p,:), [-500,  1000],  rew(p)), krnl, kw)];
    
    % spontaneous activity
    [RESP.Spont_mean(l), RESP.Spont_var(l), RESP.Spont_ste(l), ~, spontrates] = get_meanresp(VisResp(l).spiketimes, Spont);
    prct = prctile(spontrates,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.Spont_P01(l)  = prct(1);
    RESP.Spont_P02(l)  = prct(2);
    RESP.Spont_P05(l)  = prct(3);
    RESP.Spont_P10(l)  = prct(4);
    RESP.Spont_P25(l)  = prct(5);
    RESP.Spont_median(l)  = prct(6);
    RESP.Spont_P75(l)  = prct(7);
    RESP.Spont_P90(l)  = prct(8);
    RESP.Spont_P95(l)  = prct(9);
    RESP.Spont_P98(l)  = prct(10);
    RESP.Spont_P99(l)  = prct(11);
    RESP.Spont_MAD(l)  = mad(spontrates, 1);
    RESP.Spont_rate(p) = spontrates;
    
    % visual response (clipped after saccade)
    [RESP.VisTran_mean(l), RESP.VisTran_var(l), RESP.VisTran_ste(l), ~, VisTran_rate] = get_meanresp(VisResp(l).spiketimes, VisTran, sacc_corr(p)-10); % subtract 10 to account for latency
    prct = prctile(VisTran_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.VisTran_P01(l)  = prct(1);
    RESP.VisTran_P02(l)  = prct(2);
    RESP.VisTran_P05(l)  = prct(3);
    RESP.VisTran_P10(l)  = prct(4);
    RESP.VisTran_P25(l)  = prct(5);
    RESP.VisTran_median(l)  = prct(6);
    RESP.VisTran_P75(l)  = prct(7);
    RESP.VisTran_P90(l)  = prct(8);
    RESP.VisTran_P95(l)  = prct(9);
    RESP.VisTran_P98(l)  = prct(10);
    RESP.VisTran_P99(l)  = prct(11);
    RESP.VisTran_MAD(l)  = mad(VisTran_rate, 1);
    RESP.VisTran_rate(p) = VisTran_rate;
 
    [RESP.VisSust_mean(l), RESP.VisSust_var(l), RESP.VisSust_ste(l), ~, VisSust_rate] = get_meanresp(VisResp(l).spiketimes, VisSust, sacc_corr(p)-10);
    prct = prctile(VisSust_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.VisSust_P01(l)  = prct(1);
    RESP.VisSust_P02(l)  = prct(2);
    RESP.VisSust_P05(l)  = prct(3);
    RESP.VisSust_P10(l)  = prct(4);
    RESP.VisSust_P25(l)  = prct(5);
    RESP.VisSust_median(l)  = prct(6);
    RESP.VisSust_P75(l)  = prct(7);
    RESP.VisSust_P90(l)  = prct(8);
    RESP.VisSust_P95(l)  = prct(9);
    RESP.VisSust_P98(l)  = prct(10);
    RESP.VisSust_P99(l)  = prct(11);
    RESP.VisSust_MAD(l)  = mad(VisSust_rate, 1);
    RESP.VisSust_rate(p) = VisSust_rate;
 
    [RESP.VisLate_mean(l), RESP.VisLate_var(l), RESP.VisLate_ste(l), ~, VisLate_rate] = get_meanresp(VisResp(l).spiketimes, VisLate, sacc_corr(p)-10);
    prct = prctile(VisLate_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.VisLate_P01(l)  = prct(1);
    RESP.VisLate_P02(l)  = prct(2);
    RESP.VisLate_P05(l)  = prct(3);
    RESP.VisLate_P10(l)  = prct(4);
    RESP.VisLate_P25(l)  = prct(5);
    RESP.VisLate_median(l)  = prct(6);
    RESP.VisLate_P75(l)  = prct(7);
    RESP.VisLate_P90(l)  = prct(8);
    RESP.VisLate_P95(l)  = prct(9);
    RESP.VisLate_P98(l)  = prct(10);
    RESP.VisLate_P99(l)  = prct(11);
    RESP.VisLate_MAD(l)  = mad(VisLate_rate, 1);
    RESP.VisLate_rate(p) = VisLate_rate;

    % saccade related response
    [RESP.MovEarl_mean(l), RESP.MovEarl_var(l), RESP.MovEarl_ste(l), ~, MovEarl_rate] = get_meanresp(MovResp(l).spiketimes, MovEarl);
    prct = prctile(MovEarl_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.MovEarl_P01(l)  = prct(1);
    RESP.MovEarl_P02(l)  = prct(2);
    RESP.MovEarl_P05(l)  = prct(3);
    RESP.MovEarl_P10(l)  = prct(4);
    RESP.MovEarl_P25(l)  = prct(5);
    RESP.MovEarl_median(l)  = prct(6);
    RESP.MovEarl_P75(l)  = prct(7);
    RESP.MovEarl_P90(l)  = prct(8);
    RESP.MovEarl_P95(l)  = prct(9);
    RESP.MovEarl_P98(l)  = prct(10);
    RESP.MovEarl_P99(l)  = prct(11);
    RESP.MovEarl_MAD(l)  = mad(MovEarl_rate, 1);
    RESP.MovEarl_rate(p) = MovEarl_rate;
    
    [RESP.MovSust_mean(l), RESP.MovSust_var(l), RESP.MovSust_ste(l), ~, MovSust_rate] = get_meanresp(MovResp(l).spiketimes, MovSust);
    prct = prctile(MovSust_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.MovSust_P01(l)  = prct(1);
    RESP.MovSust_P02(l)  = prct(2);
    RESP.MovSust_P05(l)  = prct(3);
    RESP.MovSust_P10(l)  = prct(4);
    RESP.MovSust_P25(l)  = prct(5);
    RESP.MovSust_median(l)  = prct(6);
    RESP.MovSust_P75(l)  = prct(7);
    RESP.MovSust_P90(l)  = prct(8);
    RESP.MovSust_P95(l)  = prct(9);
    RESP.MovSust_P98(l)  = prct(10);
    RESP.MovSust_P99(l)  = prct(11);
    RESP.MovSust_MAD(l)  = mad(MovSust_rate, 1);
    RESP.MovSust_rate(p) = MovSust_rate;
    
    [RESP.MovTran_mean(l), RESP.MovTran_var(l), RESP.MovTran_ste(l), ~, MovTran_rate] = get_meanresp(MovResp(l).spiketimes, MovTran);
    prct = prctile(MovTran_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.MovTran_P01(l)  = prct(1);
    RESP.MovTran_P02(l)  = prct(2);
    RESP.MovTran_P05(l)  = prct(3);
    RESP.MovTran_P10(l)  = prct(4);
    RESP.MovTran_P25(l)  = prct(5);
    RESP.MovTran_median(l)  = prct(6);
    RESP.MovTran_P75(l)  = prct(7);
    RESP.MovTran_P90(l)  = prct(8);
    RESP.MovTran_P95(l)  = prct(9);
    RESP.MovTran_P98(l)  = prct(10);
    RESP.MovTran_P99(l)  = prct(11);
    RESP.MovTran_MAD(l)  = mad(MovTran_rate, 1);
    RESP.MovTran_rate(p) = MovTran_rate;
    
    [RESP.MovSacc_mean(l), RESP.MovSacc_var(l), RESP.MovSacc_ste(l), ~, MovSacc_rate] = get_meanresp(MovResp(l).spiketimes, MovSacc);
    prct = prctile(MovSacc_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.MovSacc_P01(l)  = prct(1);
    RESP.MovSacc_P02(l)  = prct(2);
    RESP.MovSacc_P05(l)  = prct(3);
    RESP.MovSacc_P10(l)  = prct(4);
    RESP.MovSacc_P25(l)  = prct(5);
    RESP.MovSacc_median(l)  = prct(6);
    RESP.MovSacc_P75(l)  = prct(7);
    RESP.MovSacc_P90(l)  = prct(8);
    RESP.MovSacc_P95(l)  = prct(9);
    RESP.MovSacc_P98(l)  = prct(10);
    RESP.MovSacc_P99(l)  = prct(11);
    RESP.MovSacc_MAD(l)  = mad(MovSacc_rate, 1);
    RESP.MovSacc_rate(p) = MovSacc_rate;
    
    [RESP.MovPost_mean(l), RESP.MovPost_var(l), RESP.MovPost_ste(l), ~, MovPost_rate] = get_meanresp(MovResp(l).spiketimes, MovPost);
    prct = prctile(MovPost_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.MovPost_P01(l)  = prct(1);
    RESP.MovPost_P02(l)  = prct(2);
    RESP.MovPost_P05(l)  = prct(3);
    RESP.MovPost_P10(l)  = prct(4);
    RESP.MovPost_P25(l)  = prct(5);
    RESP.MovPost_median(l)  = prct(6);
    RESP.MovPost_P75(l)  = prct(7);
    RESP.MovPost_P90(l)  = prct(8);
    RESP.MovPost_P95(l)  = prct(9);
    RESP.MovPost_P98(l)  = prct(10);
    RESP.MovPost_P99(l)  = prct(11);
    RESP.MovPost_MAD(l)  = mad(MovPost_rate, 1);
    RESP.MovPost_rate(p) = MovPost_rate;
      
    % reward related response 
    [RESP.EarlRew_mean(l), RESP.EarlRew_var(l), RESP.EarlRew_ste(l), ~, EarlRew_rate] = get_meanresp(RewResp(l).spiketimes, EarlRew); % subtract 10 to account for latency
    prct = prctile(EarlRew_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.EarlRew_P01(l)  = prct(1);
    RESP.EarlRew_P02(l)  = prct(2);
    RESP.EarlRew_P05(l)  = prct(3);
    RESP.EarlRew_P10(l)  = prct(4);
    RESP.EarlRew_P25(l)  = prct(5);
    RESP.EarlRew_median(l)  = prct(6);
    RESP.EarlRew_P75(l)  = prct(7);
    RESP.EarlRew_P90(l)  = prct(8);
    RESP.EarlRew_P95(l)  = prct(9);
    RESP.EarlRew_P98(l)  = prct(10);
    RESP.EarlRew_P99(l)  = prct(11);
    RESP.EarlRew_MAD(l)  = mad(EarlRew_rate, 1);
    RESP.EarlRew_rate(p) = EarlRew_rate;
      
    % reward related response 
    [RESP.PreRew_mean(l), RESP.PreRew_var(l), RESP.PreRew_ste(l), ~, PreRew_rate] = get_meanresp(RewResp(l).spiketimes, PreRew); % subtract 10 to account for latency
    prct = prctile(PreRew_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.PreRew_P01(l)  = prct(1);
    RESP.PreRew_P02(l)  = prct(2);
    RESP.PreRew_P05(l)  = prct(3);
    RESP.PreRew_P10(l)  = prct(4);
    RESP.PreRew_P25(l)  = prct(5);
    RESP.PreRew_median(l)  = prct(6);
    RESP.PreRew_P75(l)  = prct(7);
    RESP.PreRew_P90(l)  = prct(8);
    RESP.PreRew_P95(l)  = prct(9);
    RESP.PreRew_P98(l)  = prct(10);
    RESP.PreRew_P99(l)  = prct(11);
    RESP.PreRew_MAD(l)  = mad(PreRew_rate, 1);
    RESP.PreRew_rate(p) = PreRew_rate;
    
    [RESP.PostRew_mean(l), RESP.PostRew_var(l), RESP.PostRew_ste(l), ~, PostRew_rate] = get_meanresp(RewResp(l).spiketimes, PostRew);
    prct = prctile(PostRew_rate,[1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]);
    RESP.PostRew_P01(l)  = prct(1);
    RESP.PostRew_P02(l)  = prct(2);
    RESP.PostRew_P05(l)  = prct(3);
    RESP.PostRew_P10(l)  = prct(4);
    RESP.PostRew_P25(l)  = prct(5);
    RESP.PostRew_median(l)  = prct(6);
    RESP.PostRew_P75(l)  = prct(7);
    RESP.PostRew_P90(l)  = prct(8);
    RESP.PostRew_P95(l)  = prct(9);
    RESP.PostRew_P98(l)  = prct(10);
    RESP.PostRew_P99(l)  = prct(11);
    RESP.PostRew_MAD(l)  = mad(PostRew_rate, 1);
    RESP.PostRew_rate(p) = PostRew_rate;
    
    epoch_rates = [spontrates(:); VisTran_rate(:); VisSust_rate(:); MovEarl_rate(:); MovSust_rate(:); MovTran_rate(:); ...
                   MovSacc_rate(:); MovPost_rate(:); EarlRew_rate(:); PreRew_rate(:); PostRew_rate(:)];
               
    epoch_names = [repmat({'Spont'},length(spontrates),1); repmat({'VisTran'},length(VisTran_rate),1); ...
                   repmat({'VisSust'},length(VisSust_rate),1); repmat({'MovEarl'},length(MovEarl_rate),1); ...
                   repmat({'MovSust'},length(MovSust_rate),1); repmat({'MovSacc'},length(MovTran_rate),1); ...
                   repmat({'MovSacc'},length(MovSacc_rate),1); ...
                   repmat({'MovPost'},length(MovPost_rate),1); repmat({'EarlRew'},length(EarlRew_rate),1); ...
                   repmat({'PreRew'}, length(PreRew_rate),1);  repmat({'PostRew'},length(PostRew_rate),1)];
               
    RESP.Resp_Modulation_p = kruskalwallis(epoch_rates, epoch_names, 'off');     
    
    % ____________________________________________________________________________ %
    %% get p-values for the comparisons as well
    % use a paired Wilcoxon signed rank test, because samples are not independent
    % from each other but from the same trial.

    % visual transiency
    RESP.Vis_Tran_vs_Sust_P(l) = ranksum(VisTran_rate, VisSust_rate);
    RESP.Vis_Tran_vs_Late_P(l) = ranksum(VisTran_rate, VisLate_rate);
    RESP.Vis_Sust_vs_Late_P(l) = ranksum(VisSust_rate, VisLate_rate);

    % saccade transiency
    RESP.Mov_Tran_vs_Sust_P(l) = ranksum(MovTran_rate, MovSust_rate);
    RESP.Mov_Tran_vs_Sust_P(l) = ranksum(MovTran_rate, MovSust_rate);
    RESP.Mov_Tran_vs_Sust_P(l) = ranksum(MovTran_rate, MovSust_rate);

    RESP.Mov_Tran_vs_Sust_P(l) = ranksum(MovTran_rate, MovSust_rate);
    RESP.Mov_Tran_vs_Earl_P(l) = ranksum(MovTran_rate, MovEarl_rate);
    RESP.Mov_Sust_vs_Earl_P(l) = ranksum(MovSust_rate, MovEarl_rate);

    % saccade periods
    RESP.Mov_Sacc_vs_Post_P(l) = ranksum(MovSacc_rate, MovPost_rate);
    
    RESP.Mov_Pre_vs_Post_P(l)  = ranksum(MovTran_rate, MovPost_rate);
    RESP.Mov_Pre_vs_Sacc_P(l)  = ranksum(MovTran_rate, MovSacc_rate);
    RESP.Mov_Pre_vs_Sust_P(l)  = ranksum(MovTran_rate, MovSust_rate);
    RESP.Mov_Pre_vs_Earl_P(l)  = ranksum(MovTran_rate, MovEarl_rate);
    
    % reward periods
    RESP.RewPrePost_P(l)  = ranksum(PreRew_rate, PostRew_rate);

    % mov/vis contrast
    RESP.VisTran_vs_MovTran_P(l) = ranksum(VisTran_rate, MovTran_rate);
    RESP.VisSust_vs_MovSust_P(l) = ranksum(VisSust_rate, MovSust_rate);

    % response modulation
    RESP.VisTran_vs_Spont_P(l)   = ranksum(VisTran_rate, spontrates);
    RESP.VisSust_vs_Spont_P(l)   = ranksum(VisSust_rate, spontrates);
    RESP.MovEarl_vs_Spont_P(l)   = ranksum(MovEarl_rate, spontrates);
    RESP.MovSust_vs_Spont_P(l)   = ranksum(MovSust_rate, spontrates);
    RESP.MovTran_vs_Spont_P(l)   = ranksum(MovTran_rate, spontrates);
    RESP.MovSacc_vs_Spont_P(l)   = ranksum(MovSacc_rate, spontrates);
    RESP.PostRew_vs_Spont_P(l)   = ranksum(PostRew_rate, spontrates);    
        
%     %% get trial based modulations
%     % visual transiency
%     RESP.VisTrans_trialIDX(l)   = nanmean(get_contrast(VisTran_rate, VisSust_rate));
% 
%     % saccade transiency
%     RESP.MovTrans_trialIDX(l)   = nanmean(get_contrast(MovTran_rate, MovSust_rate));
% 
%     % saccade periods
%     RESP.MovPrePost_trialIDX(l)  = nanmean(get_contrast(MovTran_rate, MovPost_rate));
%     RESP.MovPreSacc_trialIDX(l)  = nanmean(get_contrast(MovTran_rate, MovSacc_rate));
%     RESP.MovSaccPost_trialIDX(l) = nanmean(get_contrast(MovSacc_rate, MovPost_rate));
% 
%     % reward periods
%     RESP.RewPrePost_trialIDX(l)  = nanmean(get_contrast(PreRew_rate,  PostRew_rate));
%     RESP.PostSaccRew_trialIDX(l) = nanmean(get_contrast(MovPost_rate, PostRew_rate));
% 
%     % mov/vis contrast
%     RESP.VisTran_vs_MovTran_trialIDX(l) = nanmean(get_contrast(VisTran_rate, MovTran_rate));
%     RESP.VisSust_vs_MovSust_trialIDX(l) = nanmean(get_contrast(VisSust_rate, MovSust_rate));
end

RESP.VisResp = VisResp;
RESP.MovResp = MovResp;
RESP.RewResp = RewResp;

RESP.Spont = nanmean(RESP.Spont_rate);

RESP.VisTran_norm = RESP.VisTran_mean ./ RESP.Spont;
RESP.VisSust_norm = RESP.VisSust_mean ./ RESP.Spont;

RESP.MovSust_norm = RESP.MovSust_mean ./ RESP.Spont;
RESP.MovTran_norm = RESP.MovTran_mean ./ RESP.Spont;
RESP.MovSacc_norm = RESP.MovSacc_mean ./ RESP.Spont;
RESP.MovPost_norm = RESP.MovPost_mean ./ RESP.Spont;

RESP.PreRew_norm  = RESP.PreRew_mean  ./ RESP.Spont;
RESP.PostRew_norm = RESP.PostRew_mean ./ RESP.Spont;

% ____________________________________________________________________________ %
% calculate ratios for the responses

% visual transiency
RESP.VisTrans_IDX   = get_contrast(RESP.VisTran_mean, RESP.VisSust_mean);
RESP.VisTrans2_IDX  = get_contrast(RESP.VisSust_mean, RESP.VisLate_mean);

% saccade transiency
RESP.MovTrans_IDX   = get_contrast(RESP.MovTran_mean, RESP.MovSust_mean);

% saccade periods
RESP.MovPrePost_IDX  = get_contrast(RESP.MovTran_mean, RESP.MovPost_mean);
RESP.MovPreSacc_IDX  = get_contrast(RESP.MovTran_mean, RESP.MovSacc_mean);
RESP.MovSaccPost_IDX = get_contrast(RESP.MovSacc_mean, RESP.MovPost_mean);

% reward periods
RESP.RewPrePost_IDX  = get_contrast(RESP.PreRew_mean, RESP.PostRew_mean);
RESP.RewPostSacc_IDX = get_contrast(RESP.MovPost_mean, RESP.PostRew_mean);

% mov/vis contrast
RESP.VisMov_IDX   = get_contrast(RESP.VisTran_mean, RESP.MovTran_mean);
RESP.VisSpont_IDX = get_contrast(RESP.VisTran_mean, RESP.Spont_mean);
RESP.MovSpont_IDX = get_contrast(RESP.MovTran_mean, RESP.Spont_mean);
RESP.RewSpont_IDX = get_contrast(RESP.PostRew_mean, RESP.Spont_mean);


% ____________________________________________________________________________ %
%% create plots
if(do_plot == 1)
    figure('Position', [0 0 1400 850], 'Renderer', 'Painters', 'Color', [1,1,1]);

    Vpos = VisResp(1).xtime >= -150         & VisResp(l).xtime <= sacc_prct;
    Mpos = MovResp(1).xtime >= -1*sacc_prct & MovResp(l).xtime <= 150;
    Rpos = RewResp(1).xtime >= -350         & RewResp(l).xtime <= -350+sum(Vpos); % keep the scale the same across plots

    % summary
    subplot('Position',[0.05 0.05 0.275 0.3]);
    hold on;

    med_sac = nanmedian(sacc_corr);
    med_rew = med_sac + 3*MovPost(2); % just reduce the spacing
    
    if(any(spks(:) < 0))
         shadedErrorBar(PreFix+min(Spont), RESP.PreFix_mean([1,1]), RESP.PreFix_ste([1,1]), {'color', [0.6, 0.6, 0.6],'markerfacecolor',[0.8, 0.8, 0.8]},do_alpha);
%         shadedErrorBar(Spont, RESP.PreFix_mean([1,1]), [RESP.PreFix_P05,RESP.PreFix_P95], {'color', [0.6, 0.6, 0.6],'markerfacecolor',[0.8, 0.8, 0.8]},do_alpha);
    end

    for(l=1:length(loclst))   
        % spontaneous  
        shadedErrorBar(Spont, RESP.Spont_mean([l,l]), RESP.Spont_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
%         shadedErrorBar(Spont, RESP.Spont_median([l,l]), [RESP.Spont_P05(l),RESP.Spont_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);

        % visual
        shadedErrorBar(VisTran, RESP.VisTran_mean([l,l]), RESP.VisTran_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
        shadedErrorBar(VisSust, RESP.VisSust_mean([l,l]), RESP.VisSust_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
        shadedErrorBar(VisLate, RESP.VisLate_mean([l,l]), RESP.VisLate_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
 
%         shadedErrorBar(VisTran, RESP.VisTran_median([l,l]), [RESP.VisTran_P05(l),RESP.VisTran_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
%         shadedErrorBar(VisSust, RESP.VisSust_median([l,l]), [RESP.VisSust_P05(l),RESP.VisSust_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);

        % movement
        shadedErrorBar(med_sac+MovSust, RESP.MovSust_mean([l,l]), RESP.MovSust_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},1);
        shadedErrorBar(med_sac+MovTran, RESP.MovTran_mean([l,l]), RESP.MovTran_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
        shadedErrorBar(med_sac+MovSacc, RESP.MovSacc_mean([l,l]), RESP.MovSacc_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
        shadedErrorBar(med_sac+MovPost, RESP.MovPost_mean([l,l]), RESP.MovPost_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);

%         shadedErrorBar(med_sac+MovSust, RESP.MovSust_median([l,l]), [RESP.MovSust_P05(l),RESP.MovSust_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},1);
%         shadedErrorBar(med_sac+MovTran, RESP.MovTran_median([l,l]), [RESP.MovTran_P05(l),RESP.MovTran_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
%         shadedErrorBar(med_sac+MovSacc, RESP.MovSacc_median([l,l]), [RESP.MovSacc_P05(l),RESP.MovSacc_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
%         shadedErrorBar(med_sac+MovPost, RESP.MovPost_median([l,l]), [RESP.MovPost_P05(l),RESP.MovPost_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);

        % reward
        shadedErrorBar(med_rew+EarlRew, RESP.EarlRew_mean([l,l]), RESP.EarlRew_ste([l,l]),  {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
        shadedErrorBar(med_rew+PreRew,  RESP.PreRew_mean([l,l]),  RESP.PreRew_ste([l,l]),  {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
        shadedErrorBar(med_rew+PostRew, RESP.PostRew_mean([l,l]), RESP.PostRew_ste([l,l]), {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);

%         shadedErrorBar(med_rew+EarlRew, RESP.EarlRew_median([l,l]), [RESP.EarlRew_P05(l),RESP.EarlRew_P95(l)],  {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
%         shadedErrorBar(med_rew+PreRew,  RESP.PreRew_median([l,l]),  [RESP.PreRew_P05(l),RESP.PreRew_P95(l)],  {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
%         shadedErrorBar(med_rew+PostRew, RESP.PostRew_median([l,l]), [RESP.PostRew_P05(l),RESP.PostRew_P95(l)], {'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha);
    end

    Y_lim = [0,max(ylim())];
    ylim(Y_lim);

    plot([0,0],Y_lim,'-k')
    plot([med_sac,med_sac],Y_lim,'-k')
    plot([med_rew,med_rew],Y_lim,'-k')

    text(0,1.01*Y_lim(2),'Stim',...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom','color',[0.4, 0.4,0.4]);
    text(med_sac,1.01*Y_lim(2),'Sacc',...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom','color',[0.4, 0.4,0.4]);
    text(med_rew,1.01*Y_lim(2),'Rew',...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom','color',[0.4, 0.4,0.4]);

    box off; axis tight;
    set(gca,'TickDir','out');

    % visual
    %spkplot(1) = subplot(3,3,1);
    spkplot(1) = subplot('Position',[0.05 0.75 0.275 0.2]);

    title('Visual')
    hold on;
    ctrial = 0;
    for(l=1:length(loclst))     
        p = loc == loclst(l);
        ctrial = SPK_plot_spikeraster(VisResp(l).spiketimes, sacc_corr(p), [min(VisResp(1).xtime(Vpos)), max(VisResp(1).xtime(Vpos))], ctrial, main_col(l,:), [0, 0, 0]);
    end
    xlim([min(VisResp(1).xtime(Vpos)), max(VisResp(1).xtime(Vpos))]);

    % rtplot(1) = subplot(3,3,4);
    rtplot(1) = subplot('Position',[0.05 0.45 0.275 0.3]);
    hold on;
    for(l=1:length(loclst))     
        shadedErrorBar(VisResp(l).xtime(Vpos),VisResp(l).meandensity(Vpos),VisResp(l).stedensity(Vpos),{'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha)

        if(get_rate == 1)
            if(max_rate < max(VisResp(l).meandensity(Vpos)+VisResp(l).stedensity(Vpos)))
               max_rate = max(VisResp(l).meandensity(Vpos)+VisResp(l).stedensity(Vpos));
            end
        end
    end

    % movement
    % spkplot(2) = subplot(3,3,2);
    spkplot(2) = subplot('Position',[0.375 0.75 0.275 0.2]);

    title('Movement')
    hold on;
    ctrial = 0;
    for(l=1:length(loclst))     
        p = loc == loclst(l);
        ctrial = SPK_plot_spikeraster(MovResp(l).spiketimes, -1*sacc_corr(p), [min(MovResp(1).xtime(Mpos)), max(MovResp(1).xtime(Mpos))], ctrial, main_col(l,:), [0, 0, 0]);
    end
    xlim([min(MovResp(1).xtime(Mpos)), max(MovResp(1).xtime(Mpos))]);

    % rtplot(2) = subplot(3,3,5);
    rtplot(2) = subplot('Position',[0.375 0.45 0.275 0.3]);
    hold on;
    for(l=1:length(loclst))     
        shadedErrorBar(MovResp(l).xtime(Mpos),MovResp(l).meandensity(Mpos),MovResp(l).stedensity(Mpos),{'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha)

        if(get_rate == 1)
            if(max_rate < max(MovResp(l).meandensity(Mpos)+MovResp(l).stedensity(Mpos)))
                max_rate = max(MovResp(l).meandensity(Mpos)+MovResp(l).stedensity(Mpos));
            end
        end
    end

    % reward
    % spkplot(3) = subplot(3,3,3);
    spkplot(3) = subplot('Position',[0.7 0.75 0.275 0.2]);
    title('Reward')
    hold on;
    ctrial = 0;
    for(l=1:length(loclst))     
        ctrial = SPK_plot_spikeraster(RewResp(l).spiketimes, [], [min(RewResp(1).xtime(Rpos)), max(RewResp(1).xtime(Rpos))], ctrial, main_col(l,:), [0, 0, 0]);
    end
    xlim([min(RewResp(1).xtime(Rpos)), max(RewResp(1).xtime(Rpos))]);

    % rtplot(3) = subplot(3,3,6);
    rtplot(3) = subplot('Position',[0.7 0.45 0.275 0.3]);
    hold on;
    for(l=1:length(loclst))     
        shadedErrorBar(RewResp(l).xtime(Rpos),RewResp(l).meandensity(Rpos),RewResp(l).stedensity(Rpos),{'color', main_col(l,:),'markerfacecolor',ptch_col(l,:)},do_alpha)

        if(get_rate == 1)
            if(max_rate < max(RewResp(l).meandensity(Rpos)+RewResp(l).stedensity(Rpos)))
                max_rate = max(RewResp(l).meandensity(Rpos)+RewResp(l).stedensity(Rpos));
            end
        end
    end

    % adjust scales   
    if(get_rate == 1)
        if(max_rate > 50)
            max_rate = ceil(prctile(max_rate,95)/10)*10;
        elseif(max_rate > 25)
            max_rate = ceil(prctile(max_rate,95)/5)*5;
        else
            max_rate = ceil(prctile(max_rate,95)/2)*2;
        end
    end

    if(get_trial == 1)
        max_trial = size(spks,1)+1;
    end

    for(p=1:length(rtplot))
        subplot(rtplot(p));
        box off; axis tight;
        set(gca,'TickDir','out');
        ylim([0,max_rate]);
        plot([0,0],[0,max_rate],'g', 'LineWidth',1.5);
        X_lim = xlim();

        xlabel('time [ms]')
        if(p==1)
            ylabel('firing rate [spk/s]');
        end
        
        subplot(spkplot(p));
        axis off; axis tight;
        ylim([0,max_trial]);
        xlim(X_lim);
        plot([0,0],[0,max_trial],'g', 'LineWidth',1.5);   
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions

% ____________________________________________________________________________ %
%% plot spike raster
function ctr = get_contrast(a,b)
    ctr = (a - b) ./ (a + b);

% ____________________________________________________________________________ %
%% get basic response characterisation for a defined time window
function [meanrt, varrt, stert, fano, trialrate] = get_meanresp(spktm, timwin, clip)

    windur = repmat(diff(timwin)+1,size(spktm,1),1);
    
    if(exist('clip','var') && ~isempty(clip))
        spktm(bsxfun(@gt, spktm, clip(:))) = NaN;
        
        windur2 = (clip(:) - timwin(1)) + 1;
        p = windur2 < windur;
        windur(p) = windur2(p);
    end
    
    pos = spktm >= timwin(1)  & spktm <= timwin(2);
    spkcnt    = sum(pos,2);
    trialrate = (1000 .* spkcnt) ./ windur;
    

    meanrt = nanmean(trialrate);
    varrt  =  nanvar(trialrate);
    stdrt  =  nanstd(trialrate);
    stert  =  stdrt / sqrt(sum(isfinite(trialrate)));
    fano   = varrt/meanrt;

% ____________________________________________________________________________ %
%% plot spike raster
% function ctrial = SPK_plot_spikeraster(spks, evtm, twin, startY, col, evcol)
%     
%     if(~exist('startY','var') || isempty(startY))
%         startY = 0;
%     end
%     
%     if(~exist('evtm','var'))
%         evtm = [];
%     end
%     
%     if(~exist('twin','var') || isempty(twin))
%         twin = [min(spks(:)), max(spks(:))];
%     end
%     
%     if(~exist('col','var') || isempty(col))
%         col = [0, 0, 0];
%     end
%     
%     if(~exist('evcol','var') || isempty(evcol))
%         evcol = [1, 0, 0];
%     end
%             
%     hold on;
%     
%     nTrials = size(spks,1);
%     ctrial = startY;
%     
%     if(~isempty(evtm))
%         med_evtm = nanmedian(evtm);
%         plot([med_evtm,med_evtm], [startY,startY+nTrials+1], '-', 'color', evcol, 'LineWidth', 1.5);
%         
%         [~, Torder] = sort(evtm);
%     else
%         Torder = 1:nTrials;
%     end
%     
%     if(~isempty(spks) && any(isfinite(spks(:))))
% 
%         for(s=1:nTrials)
%             ctrial  = ctrial + 1;
% 
%             cspikes = spks(Torder(s),:);
%             cspikes(isnan(cspikes)) = [];
%             cspikes(cspikes<twin(1) | cspikes>twin(2)) = [];
%             
%             if(~isempty(evtm))
%                 plot([evtm(Torder(s)),evtm(Torder(s))], [ctrial-0.5,ctrial+0.5],'color',evcol, 'LineWidth', 1.5);
%             end
%             
%             if(~isempty(cspikes))
%                 plot(cspikes,ctrial,'.', 'color', col);
%             end
% 
%         end
%     end
% 
%    