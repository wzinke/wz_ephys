function TRANS = wz_spk_getTransiency(visresp, visX, movresp, sacX, spontmat)

%  ========================================================================
%% Define constant parameter


%  ========================================================================
%% check input arguments and define defaults

% number of bootstrap repetitions


% =========================================================================
%% prepare data


%% visual transiency

TRANS.Spont = mean(spontmat(:));

TRANS.VisTran = get_meanresp(visresp, visX, [ 75; 125]);
TRANS.VisSust = get_meanresp(visresp, visX, [150; 200]);

TRANS.VisTrans_IDX_cont  = get_contrast(TRANS.VisTran, TRANS.VisSust);
TRANS.VisTrans_IDX_ratio = TRANS.VisTran / TRANS.VisSust;

%% movement transiency
TRANS.MovSust = get_meanresp(movresp, sacX, [-125; -75]);
TRANS.MovTran = get_meanresp(movresp, sacX, [ -50;   0]);
TRANS.MovSacc = get_meanresp(movresp, sacX, [   0;  50]);
TRANS.MovPost = get_meanresp(movresp, sacX, [  50; 100]);

TRANS.MovTrans_IDX_cont  = get_contrast(TRANS.MovTran, TRANS.MovSust);
TRANS.MovTrans_IDX_ratio = TRANS.MovTran / TRANS.MovSust;

TRANS.MovPrePost_IDX  = get_contrast(TRANS.MovTran, TRANS.MovPost);
TRANS.MovPreSacc_IDX  = get_contrast(TRANS.MovTran, TRANS.MovSacc);
TRANS.MovSaccPost_IDX = get_contrast(TRANS.MovSacc, TRANS.MovPost);

%% mov/vis contrast
TRANS.VisMov_IDX   = get_contrast(TRANS.VisTran, TRANS.MovTran);
TRANS.VisSpont_IDX = get_contrast(TRANS.VisTran, TRANS.Spont);
TRANS.MovSpont_IDX = get_contrast(TRANS.MovTran, TRANS.Spont);

%% norm response rates

TRANS.VisTran_norm = TRANS.VisTran / TRANS.Spont;
TRANS.VisSust_norm = TRANS.VisSust / TRANS.Spont;
TRANS.MovSust_norm = TRANS.MovSust / TRANS.Spont;
TRANS.MovTran_norm = TRANS.MovTran / TRANS.Spont;
TRANS.MovSacc_norm = TRANS.MovSacc / TRANS.Spont;
TRANS.MovPost_norm = TRANS.MovPost / TRANS.Spont;

% plot(visX,mean(visresp), 'b', 'LineWidth',2.5);
% hold on
% vline(0)
% plot(sacX,mean(movresp), 'r', 'LineWidth',2.5);
% 
% hline(mean(spontmat(:)));
% 
% xlim([-350, 350]);
% hold off;

function ctr = get_contrast(a,b)
    ctr = (a - b) / (a + b);


function meanresp = get_meanresp(rmat, xvec, timwin)
    pos = xvec >= timwin(1)  & xvec < timwin(2);
    meanresp = nanmean(nanmean(rmat(:, pos)));