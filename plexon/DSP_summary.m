function  SMRY = DSP_summary(rootdir, SessDir, odir)
% Plot spiking activity signals for all DSP data.
%
% DESCRIPTION
%
%
% SYNTAX
%
%
%   Input:
%
%       SessDir    path to session directory which has to contain an 'LFP' subdirectory
%
%       odir       Output directory to save the summary plots (default is <SessDir>/DSP/plots)
%
%
%
%   Output:
%
%       SMRY       Struct with fields for all spike descriptors
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 15-Feb-2015 by wolf zinke
%
% ToDo:  - add loop to run over several input files
%        - write summary table with basic neuron descriptors

mincorrtrials =  40; % Number of trials for a condition to be plotted
maxspkrate    =  80; % maxinum firing rate during the pre-stim period to be accepted as usefull unit.

minSRT = 60; % minimal SRT to be considered as response

ktype = 'exp';
kwdth = 20;

nBoot  = 1;
smplsz = [];

vistune = [  75 125]; % time window to determine visual tuning
movtune = [ -50   0]; % time window to determine movement tuning
sponttm = [-250   0]; % prestimulus time window for spont activity

vis_win = [-500 1000];
sac_win = [-750 750];

latwin  = [25 ; 250];

% if not specified use GUI to get the file
if(~exist('rootdir','var') || isempty(rootdir))
    rootdir = [];
elseif(ischar(SessDir))
    rootdir = [];
end

% if not specified use GUI to get the session directory
if(~exist('SessDir','var') || isempty(SessDir))
    SessDir = uigetdir('path to session data');
elseif(ischar(SessDir))
    SessDir = {SessDir}; % just keep the following code consistent in expecting a cell array
end

if(~exist('odir','var') || isempty(odir))
    do_plot = 0;
else
    do_plot = 1;
    pmkdir(odir);
end

SMRY = [];

snt = 0;
for(sess=1:length(SessDir))
    cSDir = cell2mat(SessDir(sess));

    % ____________________________________________________________________________ %
    %% get available DSP channels
    [~, cSess] = fileparts(cSDir);

    DSPlst  = dir(fullfile(rootdir, cSDir,'DSP','DSP*')); % get available DSP files

    for(d=1:length(DSPlst))

        cDSp = DSPlst(d).name;

        disp([char(cSess),'-',char(cDSp)]);

        % get available mat files
        mlst = dir(fullfile(rootdir, cSDir,'DSP',cDSp,[cSess,'_',cDSp,'_*.mat']));
        numMat = length(mlst);

        % get a list of available paradigms
        dsp_arr   = [];
        task_arr  = {};
        setsz_arr = [];
        sfx_arr   = {};
        mfl_arr   = [];

        for(j=1:numMat)
            [~,cdata] = fileparts(mlst(j).name);
            splpos    = find(cdata == '_');
            fileSFX   = cdata(splpos(2)+1:end);

            cfl = fullfile(rootdir, cSDir,'DSP',cDSp,[cSess,'_',cDSp,'_',fileSFX,'.mat']);

            dsp = load(cfl);

            % make a quick check if the firing rate is excessively high and discard if so
            spospk = sum(dsp.spiketimes(:) >= sponttm(1) & dsp.spiketimes(:) <= sponttm(2));
            spospkrt = 1000 * spospk / dsp.Ntrials /(diff(sponttm)+1);

            if(spospkrt > maxspkrate || size((~isnan(dsp.spiketimes)),2) < 4 )
                continue;
            end

            p = dsp.Task.Correct == 1;

            tsks = unique(dsp.Task.TaskType);

            for(t=1:length(tsks))
                cty = strcmp(dsp.Task.TaskType, tsks(t));
                tp = cty(:) == 1 & p(:) == 1;

                if(sum(tp) > mincorrtrials)
                    if(strcmp(tsks(t), 'MG'))
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'MG'];
                        setsz_arr = [setsz_arr; NaN];
                        sfx_arr   = [sfx_arr; fileSFX];
                        mfl_arr   = [mfl_arr; d];

                    elseif(strcmp(tsks(t), 'Cap'))

                        if(sum(dsp.Task.Singleton(tp) == 1) > mincorrtrials)
                            dsp_arr   = [dsp_arr; dsp];
                            task_arr  = [task_arr; 'Cap'];
                            setsz_arr = [setsz_arr; 0];
                            sfx_arr   = [sfx_arr; fileSFX];
                            mfl_arr   = [mfl_arr; d];
                        end

                        if(sum(dsp.Task.Singleton(tp) == 0) > mincorrtrials)
                            dsp_arr   = [dsp_arr; dsp];
                            task_arr  = [task_arr; 'Cap'];
                            setsz_arr = [setsz_arr; 1];
                            sfx_arr   = [sfx_arr; fileSFX];
                            mfl_arr   = [mfl_arr; d];
                        end

                    elseif(strcmp(tsks(t), 'Det'))
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'Det'];
                        setsz_arr = [setsz_arr; 1];
                        sfx_arr   = [sfx_arr; fileSFX];
                        mfl_arr   = [mfl_arr; d];

                    elseif(strcmp(tsks(t), 'Search'))
                        if(sum(dsp.Task.SetSize(tp) == 2) > mincorrtrials)
                            dsp_arr   = [dsp_arr; dsp];
                            task_arr  = [task_arr; 'Search'];
                            setsz_arr = [setsz_arr; 2];
                            sfx_arr   = [sfx_arr; fileSFX];
                            mfl_arr   = [mfl_arr; d];
                        end

                        if(sum(dsp.Task.SetSize(tp) == 4) > mincorrtrials)
                            dsp_arr   = [dsp_arr; dsp];
                            task_arr  = [task_arr; 'Search'];
                            setsz_arr = [setsz_arr; 4];
                            sfx_arr   = [sfx_arr; fileSFX];
                            mfl_arr   = [mfl_arr; d];
                        end

                        if(sum(dsp.Task.SetSize(tp) == 8) > mincorrtrials)
                            dsp_arr   = [dsp_arr; dsp];
                            task_arr  = [task_arr; 'Search'];
                            setsz_arr = [setsz_arr; 8];
                            sfx_arr   = [sfx_arr; fileSFX];
                            mfl_arr   = [mfl_arr; d];
                        end
                    end
                end
            end  %   for(t=1:length(tsks))
        end  %  for(j=1:numMat)

        if(isempty(dsp_arr))
            continue;
        end

        % ____________________________________________________________________________ %
        %% create plot outline:  Wave | SRT Distri | Stim on | Sacc on | Reward | Error | LvsR Stim | LvsR Sacc | Ori Stim | Ori Sacc
        num_rows = length(dsp_arr);
        num_cols = 10;

        figsz  = [0 0 num_cols*360 num_rows*200];
        fhndl  = figure('Name', [cSess, '  -  ', cDSp], 'Position', figsz, 'Renderer', 'Painters');

        Ystart = linspace(98.5,4,num_rows+1)./100;
        Ystart(1)   = [];

        Xstart = linspace(4,98.5,num_cols+1)./100;
        Xstart(end) = [];

        xwd = min(diff(Xstart))-0.025;
        if(num_rows>1)
            ywd = min(abs(diff(Ystart)))-0.05;
        else
            ywd = 0.9;
        end

        plthdl = nan(num_rows,num_cols);

        for(y=1:num_rows)
            cY = Ystart(y);
            for(x=1:num_cols)
                cX = Xstart(x);

                plthdl(y,x) = subplot('Position',[cX cY xwd ywd]);
                hold on; box on; axis tight;
                set(gca,'TickDir','out');
            end
        end

        % loop over the files for each DSP channel within a session
        for(j=1:num_rows)

            disp(['...', char(sfx_arr(j)), ' - ', char(task_arr(j))]);

            spk = dsp_arr(j);
            snt = snt+1;

            %% get main info for this unit and file
            SMRY.NHP(snt,1)          = spk.NHP;
            SMRY.Date(snt,1)         = {spk.Date};
            SMRY.Paradigm(snt,1)     = {spk.Paradigm};
            SMRY.Spec(snt,1)         = {spk.Spec};
            SMRY.sessID(snt,1)       = {spk.sessID};
            SMRY.unitID(snt,1)       = {[spk.sessID,'_',spk.DSPname]};
            SMRY.channel(snt,1)      = spk.channel;
            SMRY.area(snt,1)         = {spk.area};
            SMRY.ml_coor(snt,1)      = {spk.ml_coor};
            SMRY.ap_coor(snt,1)      = {spk.ap_coor};
            SMRY.depth_uncorr(snt,1) = spk.depth_uncorr;
            SMRY.DSPname(snt,1)      = {spk.DSPname};
            SMRY.Task(snt,1)         = task_arr(j);
            SMRY.FileNum(snt,1)      = mfl_arr(j);
            SMRY.File(snt,1)         = {spk.File};

           %% get relevant trial positions

            tskp = strcmp(spk.Task.TaskType, task_arr(j));
            srtp = spk.Task.SRT > minSRT;
            szp = ones(size(tskp));

            switch char(task_arr(j))
                case 'MG'
                    detstr = [];
                    SMRY.SetSize(snt,1)   = 1;
                    SMRY.Singleton(snt,1) = 0;

                case {'Search', 'Det'}
                    szp = spk.Task.SetSize ==  setsz_arr(j) & spk.Task.IsCatch == 0;
                    if(setsz_arr(j) == 1)
                        detstr = [' - Detection'];
                    else
                        detstr = [' - Set Size ',num2str(setsz_arr(j))];
                    end
                    SMRY.SetSize(snt,1)   = setsz_arr(j);
                    SMRY.Singleton(snt,1) = 0;

                case 'Cap'
                    szp = spk.Task.Singleton == setsz_arr(j) & spk.Task.IsCatch == 0;
                    if(setsz_arr(j) == 1)
                        detstr = [' - With Singleton'];
                    else
                        detstr = [' - No Singleton'];
                    end
                    SMRY.SetSize(snt,1)   = 8;
                    SMRY.Singleton(snt,1) = setsz_arr(j);
            end

            ccnd  = tskp(:) == 1 & srtp(:) == 1 & szp(:) == 1;
            hit   = find(ccnd == 1 & spk.Task.Correct == 1); % select correct trials
            false = find(ccnd == 1 & spk.Task.error   == 1); % select false response
            left  = spk.Task.TargetLoc(hit) == 180;
            right = spk.Task.TargetLoc(hit) ==   0;

            SMRY.NumTargetLoc(snt,1) = length(unique(spk.Task.TargetLoc(hit)));

            %% get spike objects
            align_Clip = wz_spk_density(wz_get_SPKobj(spk.spiketimes(hit,:),   vis_win, [], spk.Task.SRT(hit)),    ktype, kwdth, [], spk.Task.SRT(hit));
            align_Sacc = wz_spk_density(wz_get_SPKobj(spk.spiketimes(hit,:),   sac_win,     spk.Task.SRT(hit)),    ktype, kwdth);
            align_Rew  = wz_spk_density(wz_get_SPKobj(spk.spiketimes(hit,:),   sac_win,     spk.Task.Reward(hit)), ktype, kwdth);
            align_fFB  = wz_spk_density(wz_get_SPKobj(spk.spiketimes(false,:), sac_win,     spk.Task.ErrorTone(false)), ktype, kwdth);

            vis_plw = [-100 ; prctile(spk.Task.SRT(hit),90)];
            sac_plw = [-1 *   prctile(spk.Task.SRT(hit),90); 100];

            stimpl = align_Clip.xtime > vis_plw(1) & align_Clip.xtime < vis_plw(2);
            sacmpl = align_fFB.xtime  > sac_plw(1) & align_fFB.xtime  < sac_plw(2);

            if(~isempty(align_fFB.upperdense))
                maxrate = 1.1*max( [nanmax(align_Clip.upperdense(stimpl)); ...
                    nanmax(align_Sacc.upperdense(sacmpl)); ...
                    nanmax( align_Rew.upperdense(sacmpl)); ...
                    nanmax( align_fFB.upperdense(sacmpl))]);
            else
                maxrate = 1.1*max( [nanmax(align_Clip.upperdense(stimpl)); ...
                    nanmax(align_Sacc.upperdense(sacmpl)); ...
                    nanmax( align_Rew.upperdense(sacmpl))]);
            end

            if(j==1)
                nfo = [cSess, '  -  ', cDSp];
            else
                nfo = [];
            end

            % split response across target location
            spntwin = align_Clip.xtime > sponttm(1) & align_Clip.xtime < sponttm(2);
            spontrrate = mean(nanmean(align_Clip.spikedensities(:,spntwin)));
            spontrerr  =  std(nanmean(align_Clip.spikedensities(:,spntwin),2));

            trgtloc = sort(unique(spk.Task.TargetLoc(hit)));
            thv     = deg2rad(trgtloc');

            stimtm = align_Clip.xtime > vistune(1) & align_Clip.xtime < vistune(2);
            sacctm = align_Sacc.xtime > movtune(1) & align_Sacc.xtime < movtune(2);

            mean_stim = nan(1,length(trgtloc));
            peak_stim = nan(1,length(trgtloc));
            mean_sacc = nan(1,length(trgtloc));
            peak_sacc = nan(1,length(trgtloc));

            for(a=1:length(trgtloc));
                tlpos = spk.Task.TargetLoc(hit) == trgtloc(a);

                mean_stim(a) = mean(nanmean(align_Clip.spikedensities(tlpos,stimtm)));
                peak_stim(a) =  max(nanmean(align_Clip.spikedensities(tlpos,stimtm)));

                mean_sacc(a) = mean(nanmean(align_Sacc.spikedensities(tlpos,sacctm)));
                peak_sacc(a) =  max(nanmean(align_Sacc.spikedensities(tlpos,sacctm)));
            end
            
            %% get preferred location
            prefpos_vis = trgtloc(find(mean_stim == max(mean_stim),1,'first'));
            if(~isempty(prefpos_vis))
                SMRY.PrefPos_Vis(snt,1) = prefpos_vis;
            else
                SMRY.PrefPos_Vis(snt,1) = NaN;
            end
            prefvis = spk.Task.TargetLoc(hit) == prefpos_vis;
            
            prefpos_mov = trgtloc(find(mean_sacc == max(mean_sacc),1,'first'));
            if(~isempty(prefpos_mov))
                SMRY.PrefPos_Mov(snt,1) = prefpos_mov;
            else
                SMRY.PrefPos_Mov(snt,1) = NaN;
            end
            prefmov = spk.Task.TargetLoc(hit) == prefpos_mov;


            %% get performance numbers into summary struct
            SRT_hit   = spk.Task.SRT(hit);
            SRT_false = spk.Task.SRT(false);
            SRT_left  = SRT_hit(left);
            SRT_right = SRT_hit(right);
            SRT_prefvis = SRT_hit(prefvis);
            SRT_prefmov = SRT_hit(prefmov);

            SMRY.Ntrials(snt,1)    = sum(ccnd);
            SMRY.Nhit(snt,1)       = length(hit);
            SMRY.Nfalse(snt,1)     = length(false);
            SMRY.Nright(snt,1)     = sum(left);
            SMRY.Nleft(snt,1)      = sum(right);

            SMRY.medSRT_hit(snt,1) = median(SRT_hit);
            SMRY.medSRT_(snt,1)    = median(SRT_false);
            SMRY.medSRT_(snt,1)    = median(SRT_left);
            SMRY.medSRT_(snt,1)    = median(SRT_right);

           %% response profile quantification
            RESP_hit   = SPK_resp_profile(spk.spiketimes(hit,:), SRT_hit,1);
            SMRY.Spont_mean_hit(snt,1)       = RESP_hit.Spont_mean;
            SMRY.Spont_var_hit(snt,1)        = RESP_hit.Spont_var;
            SMRY.Spont_fano_hit(snt,1)       = RESP_hit.Spont_fano;
            SMRY.VisTran_mean_hit(snt,1)     = RESP_hit.VisTran_mean;
            SMRY.VisTran_var_hit(snt,1)      = RESP_hit.VisTran_var;
            SMRY.VisTran_fano_hit(snt,1)     = RESP_hit.VisTran_fano;
            SMRY.VisSust_mean_hit(snt,1)     = RESP_hit.VisSust_mean;
            SMRY.VisSust_var_hit(snt,1)      = RESP_hit.VisSust_var;
            SMRY.VisSust_fano_hit(snt,1)     = RESP_hit.VisSust_fano;
            SMRY.MovSust_mean_hit(snt,1)     = RESP_hit.MovSust_mean;
            SMRY.MovSust_var_hit(snt,1)      = RESP_hit.MovSust_var;
            SMRY.MovSust_fano_hit(snt,1)     = RESP_hit.MovSust_fano;
            SMRY.MovTran_mean_hit(snt,1)     = RESP_hit.MovTran_mean;
            SMRY.MovTran_var_hit(snt,1)      = RESP_hit.MovTran_var;
            SMRY.MovTran_fano_hit(snt,1)     = RESP_hit.MovTran_fano;
            SMRY.MovSacc_mean_hit(snt,1)     = RESP_hit.MovSacc_mean;
            SMRY.MovSacc_var_hit(snt,1)      = RESP_hit.MovSacc_var;
            SMRY.MovSacc_fano_hit(snt,1)     = RESP_hit.MovSacc_fano;
            SMRY.MovPost_mean_hit(snt,1)     = RESP_hit.MovPost_mean;
            SMRY.MovPost_var_hit(snt,1)      = RESP_hit.MovPost_var;
            SMRY.MovPost_fano_hit(snt,1)     = RESP_hit.MovPost_fano;
            SMRY.VisResp_hit(snt,1)          = RESP_hit.VisResp;
            SMRY.MovResp_hit(snt,1)          = RESP_hit.MovResp;
            SMRY.resptype_hit(snt,1)         = {RESP_hit.resptype};
            SMRY.VisTrans_IDX_hit(snt,1)     = RESP_hit.VisTrans_IDX;
            SMRY.MovTrans_IDX_hit(snt,1)     = RESP_hit.MovTrans_IDX;
            SMRY.MovPrePost_IDX_hit(snt,1)   = RESP_hit.MovPrePost_IDX;
            SMRY.MovPreSacc_IDX_hit(snt,1)   = RESP_hit.MovPreSacc_IDX;
            SMRY.MovSaccPost_IDX_hit(snt,1)  = RESP_hit.MovSaccPost_IDX;
            SMRY.VisMov_IDX_hit(snt,1)       = RESP_hit.VisMov_IDX;
            SMRY.VisSpont_IDX_hit(snt,1)     = RESP_hit.VisSpont_IDX;
            SMRY.MovSpont_IDX_hit(snt,1)     = RESP_hit.MovSpont_IDX;
            SMRY.VisTrans_P_hit(snt,1)       = RESP_hit.VisTrans_P;
            SMRY.MovTrans_P_hit(snt,1)       = RESP_hit.MovTrans_P;
            SMRY.MovPrePost_P_hit(snt,1)     = RESP_hit.MovPrePost_P;
            SMRY.MovPreSacc_P_hit(snt,1)     = RESP_hit.MovPreSacc_P;
            SMRY.MovSaccPost_P_hit(snt,1)    = RESP_hit.MovSaccPost_P;
            SMRY.VisMov_P_hit(snt,1)         = RESP_hit.VisMov_P;
            SMRY.VisSpont_P_hit(snt,1)       = RESP_hit.VisSpont_P;
            SMRY.MovSpont_P_hit(snt,1)       = RESP_hit.MovSpont_P;
            SMRY.VisTran_Mod_hit(snt,1)      = RESP_hit.VisTran_Mod;
            SMRY.VisSust_Mod_hit(snt,1)      = RESP_hit.VisSust_Mod;
            SMRY.MovTran_Mod_hit(snt,1)      = RESP_hit.MovTran_Mod;
            SMRY.MovSacc_Mod_hit(snt,1)      = RESP_hit.MovSacc_Mod;
            SMRY.MovPost_Mod_hit(snt,1)      = RESP_hit.MovPost_Mod;

            RESP_left  = SPK_resp_profile(spk.spiketimes(hit(left),:), SRT_left,1);
            SMRY.Spont_mean_left(snt,1)       = RESP_left.Spont_mean;
            SMRY.Spont_var_left(snt,1)        = RESP_left.Spont_var;
            SMRY.Spont_fano_left(snt,1)       = RESP_left.Spont_fano;
            SMRY.VisTran_mean_left(snt,1)     = RESP_left.VisTran_mean;
            SMRY.VisTran_var_left(snt,1)      = RESP_left.VisTran_var;
            SMRY.VisTran_fano_left(snt,1)     = RESP_left.VisTran_fano;
            SMRY.VisSust_mean_left(snt,1)     = RESP_left.VisSust_mean;
            SMRY.VisSust_var_left(snt,1)      = RESP_left.VisSust_var;
            SMRY.VisSust_fano_left(snt,1)     = RESP_left.VisSust_fano;
            SMRY.MovSust_mean_left(snt,1)     = RESP_left.MovSust_mean;
            SMRY.MovSust_var_left(snt,1)      = RESP_left.MovSust_var;
            SMRY.MovSust_fano_left(snt,1)     = RESP_left.MovSust_fano;
            SMRY.MovTran_mean_left(snt,1)     = RESP_left.MovTran_mean;
            SMRY.MovTran_var_left(snt,1)      = RESP_left.MovTran_var;
            SMRY.MovTran_fano_left(snt,1)     = RESP_left.MovTran_fano;
            SMRY.MovSacc_mean_left(snt,1)     = RESP_left.MovSacc_mean;
            SMRY.MovSacc_var_left(snt,1)      = RESP_left.MovSacc_var;
            SMRY.MovSacc_fano_left(snt,1)     = RESP_left.MovSacc_fano;
            SMRY.MovPost_mean_left(snt,1)     = RESP_left.MovPost_mean;
            SMRY.MovPost_var_left(snt,1)      = RESP_left.MovPost_var;
            SMRY.MovPost_fano_left(snt,1)     = RESP_left.MovPost_fano;
            SMRY.VisResp_left(snt,1)          = RESP_left.VisResp;
            SMRY.MovResp_left(snt,1)          = RESP_left.MovResp;
            SMRY.resptype_left(snt,1)         = {RESP_left.resptype};
            SMRY.VisTrans_IDX_left(snt,1)     = RESP_left.VisTrans_IDX;
            SMRY.MovTrans_IDX_left(snt,1)     = RESP_left.MovTrans_IDX;
            SMRY.MovPrePost_IDX_left(snt,1)   = RESP_left.MovPrePost_IDX;
            SMRY.MovPreSacc_IDX_left(snt,1)   = RESP_left.MovPreSacc_IDX;
            SMRY.MovSaccPost_IDX_left(snt,1)  = RESP_left.MovSaccPost_IDX;
            SMRY.VisMov_IDX_left(snt,1)       = RESP_left.VisMov_IDX;
            SMRY.VisSpont_IDX_left(snt,1)     = RESP_left.VisSpont_IDX;
            SMRY.MovSpont_IDX_left(snt,1)     = RESP_left.MovSpont_IDX;
            SMRY.VisTrans_P_left(snt,1)       = RESP_left.VisTrans_P;
            SMRY.MovTrans_P_left(snt,1)       = RESP_left.MovTrans_P;
            SMRY.MovPrePost_P_left(snt,1)     = RESP_left.MovPrePost_P;
            SMRY.MovPreSacc_P_left(snt,1)     = RESP_left.MovPreSacc_P;
            SMRY.MovSaccPost_P_left(snt,1)    = RESP_left.MovSaccPost_P;
            SMRY.VisMov_P_left(snt,1)         = RESP_left.VisMov_P;
            SMRY.VisSpont_P_left(snt,1)       = RESP_left.VisSpont_P;
            SMRY.MovSpont_P_left(snt,1)       = RESP_left.MovSpont_P;
            SMRY.VisTran_Mod_left(snt,1)      = RESP_left.VisTran_Mod;
            SMRY.VisSust_Mod_left(snt,1)      = RESP_left.VisSust_Mod;
            SMRY.MovTran_Mod_left(snt,1)      = RESP_left.MovTran_Mod;
            SMRY.MovSacc_Mod_left(snt,1)      = RESP_left.MovSacc_Mod;
            SMRY.MovPost_Mod_left(snt,1)      = RESP_left.MovPost_Mod;

            RESP_right = SPK_resp_profile(spk.spiketimes(hit(right),:), SRT_right,1);
            SMRY.Spont_mean_right(snt,1)       = RESP_right.Spont_mean;
            SMRY.Spont_var_right(snt,1)        = RESP_right.Spont_var;
            SMRY.Spont_fano_right(snt,1)       = RESP_right.Spont_fano;
            SMRY.VisTran_mean_right(snt,1)     = RESP_right.VisTran_mean;
            SMRY.VisTran_var_right(snt,1)      = RESP_right.VisTran_var;
            SMRY.VisTran_fano_right(snt,1)     = RESP_right.VisTran_fano;
            SMRY.VisSust_mean_right(snt,1)     = RESP_right.VisSust_mean;
            SMRY.VisSust_var_right(snt,1)      = RESP_right.VisSust_var;
            SMRY.VisSust_fano_right(snt,1)     = RESP_right.VisSust_fano;
            SMRY.MovSust_mean_right(snt,1)     = RESP_right.MovSust_mean;
            SMRY.MovSust_var_right(snt,1)      = RESP_right.MovSust_var;
            SMRY.MovSust_fano_right(snt,1)     = RESP_right.MovSust_fano;
            SMRY.MovTran_mean_right(snt,1)     = RESP_right.MovTran_mean;
            SMRY.MovTran_var_right(snt,1)      = RESP_right.MovTran_var;
            SMRY.MovTran_fano_right(snt,1)     = RESP_right.MovTran_fano;
            SMRY.MovSacc_mean_right(snt,1)     = RESP_right.MovSacc_mean;
            SMRY.MovSacc_var_right(snt,1)      = RESP_right.MovSacc_var;
            SMRY.MovSacc_fano_right(snt,1)     = RESP_right.MovSacc_fano;
            SMRY.MovPost_mean_right(snt,1)     = RESP_right.MovPost_mean;
            SMRY.MovPost_var_right(snt,1)      = RESP_right.MovPost_var;
            SMRY.MovPost_fano_right(snt,1)     = RESP_right.MovPost_fano;
            SMRY.VisResp_right(snt,1)          = RESP_right.VisResp;
            SMRY.MovResp_right(snt,1)          = RESP_right.MovResp;
            SMRY.resptype_right(snt,1)         = {RESP_right.resptype};
            SMRY.VisTrans_IDX_right(snt,1)     = RESP_right.VisTrans_IDX;
            SMRY.MovTrans_IDX_right(snt,1)     = RESP_right.MovTrans_IDX;
            SMRY.MovPrePost_IDX_right(snt,1)   = RESP_right.MovPrePost_IDX;
            SMRY.MovPreSacc_IDX_right(snt,1)   = RESP_right.MovPreSacc_IDX;
            SMRY.MovSaccPost_IDX_right(snt,1)  = RESP_right.MovSaccPost_IDX;
            SMRY.VisMov_IDX_right(snt,1)       = RESP_right.VisMov_IDX;
            SMRY.VisSpont_IDX_right(snt,1)     = RESP_right.VisSpont_IDX;
            SMRY.MovSpont_IDX_right(snt,1)     = RESP_right.MovSpont_IDX;
            SMRY.VisTrans_P_right(snt,1)       = RESP_right.VisTrans_P;
            SMRY.MovTrans_P_right(snt,1)       = RESP_right.MovTrans_P;
            SMRY.MovPrePost_P_right(snt,1)     = RESP_right.MovPrePost_P;
            SMRY.MovPreSacc_P_right(snt,1)     = RESP_right.MovPreSacc_P;
            SMRY.MovSaccPost_P_right(snt,1)    = RESP_right.MovSaccPost_P;
            SMRY.VisMov_P_right(snt,1)         = RESP_right.VisMov_P;
            SMRY.VisSpont_P_right(snt,1)       = RESP_right.VisSpont_P;
            SMRY.MovSpont_P_right(snt,1)       = RESP_right.MovSpont_P;
            SMRY.VisTran_Mod_right(snt,1)      = RESP_right.VisTran_Mod;
            SMRY.VisSust_Mod_right(snt,1)      = RESP_right.VisSust_Mod;
            SMRY.MovTran_Mod_right(snt,1)      = RESP_right.MovTran_Mod;
            SMRY.MovSacc_Mod_right(snt,1)      = RESP_right.MovSacc_Mod;
            SMRY.MovPost_Mod_right(snt,1)      = RESP_right.MovPost_Mod;
            
            RESP_PrefVis = SPK_resp_profile(spk.spiketimes(hit(prefvis),:), SRT_prefvis,1);
            SMRY.Spont_mean_PrefVis(snt,1)       = RESP_PrefVis.Spont_mean;
            SMRY.Spont_var_PrefVis(snt,1)        = RESP_PrefVis.Spont_var;
            SMRY.Spont_fano_PrefVis(snt,1)       = RESP_PrefVis.Spont_fano;
            SMRY.VisTran_mean_PrefVis(snt,1)     = RESP_PrefVis.VisTran_mean;
            SMRY.VisTran_var_PrefVis(snt,1)      = RESP_PrefVis.VisTran_var;
            SMRY.VisTran_fano_PrefVis(snt,1)     = RESP_PrefVis.VisTran_fano;
            SMRY.VisSust_mean_PrefVis(snt,1)     = RESP_PrefVis.VisSust_mean;
            SMRY.VisSust_var_PrefVis(snt,1)      = RESP_PrefVis.VisSust_var;
            SMRY.VisSust_fano_PrefVis(snt,1)     = RESP_PrefVis.VisSust_fano;
            SMRY.MovSust_mean_PrefVis(snt,1)     = RESP_PrefVis.MovSust_mean;
            SMRY.MovSust_var_PrefVis(snt,1)      = RESP_PrefVis.MovSust_var;
            SMRY.MovSust_fano_PrefVis(snt,1)     = RESP_PrefVis.MovSust_fano;
            SMRY.MovTran_mean_PrefVis(snt,1)     = RESP_PrefVis.MovTran_mean;
            SMRY.MovTran_var_PrefVis(snt,1)      = RESP_PrefVis.MovTran_var;
            SMRY.MovTran_fano_PrefVis(snt,1)     = RESP_PrefVis.MovTran_fano;
            SMRY.MovSacc_mean_PrefVis(snt,1)     = RESP_PrefVis.MovSacc_mean;
            SMRY.MovSacc_var_PrefVis(snt,1)      = RESP_PrefVis.MovSacc_var;
            SMRY.MovSacc_fano_PrefVis(snt,1)     = RESP_PrefVis.MovSacc_fano;
            SMRY.MovPost_mean_PrefVis(snt,1)     = RESP_PrefVis.MovPost_mean;
            SMRY.MovPost_var_PrefVis(snt,1)      = RESP_PrefVis.MovPost_var;
            SMRY.MovPost_fano_PrefVis(snt,1)     = RESP_PrefVis.MovPost_fano;
            SMRY.VisResp_PrefVis(snt,1)          = RESP_PrefVis.VisResp;
            SMRY.MovResp_PrefVis(snt,1)          = RESP_PrefVis.MovResp;
            SMRY.resptype_PrefVis(snt,1)         = {RESP_PrefVis.resptype};
            SMRY.VisTrans_IDX_PrefVis(snt,1)     = RESP_PrefVis.VisTrans_IDX;
            SMRY.MovTrans_IDX_PrefVis(snt,1)     = RESP_PrefVis.MovTrans_IDX;
            SMRY.MovPrePost_IDX_PrefVis(snt,1)   = RESP_PrefVis.MovPrePost_IDX;
            SMRY.MovPreSacc_IDX_PrefVis(snt,1)   = RESP_PrefVis.MovPreSacc_IDX;
            SMRY.MovSaccPost_IDX_PrefVis(snt,1)  = RESP_PrefVis.MovSaccPost_IDX;
            SMRY.VisMov_IDX_PrefVis(snt,1)       = RESP_PrefVis.VisMov_IDX;
            SMRY.VisSpont_IDX_PrefVis(snt,1)     = RESP_PrefVis.VisSpont_IDX;
            SMRY.MovSpont_IDX_PrefVis(snt,1)     = RESP_PrefVis.MovSpont_IDX;
            SMRY.VisTrans_P_PrefVis(snt,1)       = RESP_PrefVis.VisTrans_P;
            SMRY.MovTrans_P_PrefVis(snt,1)       = RESP_PrefVis.MovTrans_P;
            SMRY.MovPrePost_P_PrefVis(snt,1)     = RESP_PrefVis.MovPrePost_P;
            SMRY.MovPreSacc_P_PrefVis(snt,1)     = RESP_PrefVis.MovPreSacc_P;
            SMRY.MovSaccPost_P_PrefVis(snt,1)    = RESP_PrefVis.MovSaccPost_P;
            SMRY.VisMov_P_PrefVis(snt,1)         = RESP_PrefVis.VisMov_P;
            SMRY.VisSpont_P_PrefVis(snt,1)       = RESP_PrefVis.VisSpont_P;
            SMRY.MovSpont_P_PrefVis(snt,1)       = RESP_PrefVis.MovSpont_P;
            SMRY.VisTran_Mod_PrefVis(snt,1)      = RESP_PrefVis.VisTran_Mod;
            SMRY.VisSust_Mod_PrefVis(snt,1)      = RESP_PrefVis.VisSust_Mod;
            SMRY.MovTran_Mod_PrefVis(snt,1)      = RESP_PrefVis.MovTran_Mod;
            SMRY.MovSacc_Mod_PrefVis(snt,1)      = RESP_PrefVis.MovSacc_Mod;
            SMRY.MovPost_Mod_PrefVis(snt,1)      = RESP_PrefVis.MovPost_Mod;          
            
            RESP_PrefMov = SPK_resp_profile(spk.spiketimes(hit(prefmov),:), SRT_prefmov,1);
            SMRY.Spont_mean_PrefMov(snt,1)       = RESP_PrefMov.Spont_mean;
            SMRY.Spont_var_PrefMov(snt,1)        = RESP_PrefMov.Spont_var;
            SMRY.Spont_fano_PrefMov(snt,1)       = RESP_PrefMov.Spont_fano;
            SMRY.VisTran_mean_PrefMov(snt,1)     = RESP_PrefMov.VisTran_mean;
            SMRY.VisTran_var_PrefMov(snt,1)      = RESP_PrefMov.VisTran_var;
            SMRY.VisTran_fano_PrefMov(snt,1)     = RESP_PrefMov.VisTran_fano;
            SMRY.VisSust_mean_PrefMov(snt,1)     = RESP_PrefMov.VisSust_mean;
            SMRY.VisSust_var_PrefMov(snt,1)      = RESP_PrefMov.VisSust_var;
            SMRY.VisSust_fano_PrefMov(snt,1)     = RESP_PrefMov.VisSust_fano;
            SMRY.MovSust_mean_PrefMov(snt,1)     = RESP_PrefMov.MovSust_mean;
            SMRY.MovSust_var_PrefMov(snt,1)      = RESP_PrefMov.MovSust_var;
            SMRY.MovSust_fano_PrefMov(snt,1)     = RESP_PrefMov.MovSust_fano;
            SMRY.MovTran_mean_PrefMov(snt,1)     = RESP_PrefMov.MovTran_mean;
            SMRY.MovTran_var_PrefMov(snt,1)      = RESP_PrefMov.MovTran_var;
            SMRY.MovTran_fano_PrefMov(snt,1)     = RESP_PrefMov.MovTran_fano;
            SMRY.MovSacc_mean_PrefMov(snt,1)     = RESP_PrefMov.MovSacc_mean;
            SMRY.MovSacc_var_PrefMov(snt,1)      = RESP_PrefMov.MovSacc_var;
            SMRY.MovSacc_fano_PrefMov(snt,1)     = RESP_PrefMov.MovSacc_fano;
            SMRY.MovPost_mean_PrefMov(snt,1)     = RESP_PrefMov.MovPost_mean;
            SMRY.MovPost_var_PrefMov(snt,1)      = RESP_PrefMov.MovPost_var;
            SMRY.MovPost_fano_PrefMov(snt,1)     = RESP_PrefMov.MovPost_fano;
            SMRY.VisResp_PrefMov(snt,1)          = RESP_PrefMov.VisResp;
            SMRY.MovResp_PrefMov(snt,1)          = RESP_PrefMov.MovResp;
            SMRY.resptype_PrefMov(snt,1)         = {RESP_PrefMov.resptype};
            SMRY.VisTrans_IDX_PrefMov(snt,1)     = RESP_PrefMov.VisTrans_IDX;
            SMRY.MovTrans_IDX_PrefMov(snt,1)     = RESP_PrefMov.MovTrans_IDX;
            SMRY.MovPrePost_IDX_PrefMov(snt,1)   = RESP_PrefMov.MovPrePost_IDX;
            SMRY.MovPreSacc_IDX_PrefMov(snt,1)   = RESP_PrefMov.MovPreSacc_IDX;
            SMRY.MovSaccPost_IDX_PrefMov(snt,1)  = RESP_PrefMov.MovSaccPost_IDX;
            SMRY.VisMov_IDX_PrefMov(snt,1)       = RESP_PrefMov.VisMov_IDX;
            SMRY.VisSpont_IDX_PrefMov(snt,1)     = RESP_PrefMov.VisSpont_IDX;
            SMRY.MovSpont_IDX_PrefMov(snt,1)     = RESP_PrefMov.MovSpont_IDX;
            SMRY.VisTrans_P_PrefMov(snt,1)       = RESP_PrefMov.VisTrans_P;
            SMRY.MovTrans_P_PrefMov(snt,1)       = RESP_PrefMov.MovTrans_P;
            SMRY.MovPrePost_P_PrefMov(snt,1)     = RESP_PrefMov.MovPrePost_P;
            SMRY.MovPreSacc_P_PrefMov(snt,1)     = RESP_PrefMov.MovPreSacc_P;
            SMRY.MovSaccPost_P_PrefMov(snt,1)    = RESP_PrefMov.MovSaccPost_P;
            SMRY.VisMov_P_PrefMov(snt,1)         = RESP_PrefMov.VisMov_P;
            SMRY.VisSpont_P_PrefMov(snt,1)       = RESP_PrefMov.VisSpont_P;
            SMRY.MovSpont_P_PrefMov(snt,1)       = RESP_PrefMov.MovSpont_P;
            SMRY.VisTran_Mod_PrefMov(snt,1)      = RESP_PrefMov.VisTran_Mod;
            SMRY.VisSust_Mod_PrefMov(snt,1)      = RESP_PrefMov.VisSust_Mod;
            SMRY.MovTran_Mod_PrefMov(snt,1)      = RESP_PrefMov.MovTran_Mod;
            SMRY.MovSacc_Mod_PrefMov(snt,1)      = RESP_PrefMov.MovSacc_Mod;
            SMRY.MovPost_Mod_PrefMov(snt,1)      = RESP_PrefMov.MovPost_Mod;
                        
            %% response latency
            try
                [~, latObj_hit] = wz_spk_latency(spk.spiketimes(hit,:), latwin, sponttm, nBoot, smplsz,[],[],0);

                SMRY.LatFirstPeak_hit(snt,1) =  latObj_hit.first_peak;
                SMRY.LatMaxPeak_hit(snt,1)   =  latObj_hit.max_peak;
                SMRY.LatFirstResp_hit(snt,1) =  latObj_hit.peakresp(1);
                SMRY.LatMaxResp_hit(snt,1)   =  latObj_hit.peakresp(2);
                SMRY.LatHalfMax_hit(snt,1)   =  latObj_hit.halfmax;
                SMRY.LatMaxPeak_hit(snt,1)   =  latObj_hit.maxpeak;
                SMRY.LatFirstPeak_hit(snt,1) =  latObj_hit.frstpeak;
                SMRY.LatPEP_hit(snt,1)       =  latObj_hit.PEP;
                SMRY.LatTrise_hit(snt,1)     =  latObj_hit.Trise;
            catch
                SMRY.LatFirstPeak_hit(snt,1) =  NaN;
                SMRY.LatMaxPeak_hit(snt,1)   =  NaN;
                SMRY.LatFirstResp_hit(snt,1) =  NaN;
                SMRY.LatMaxResp_hit(snt,1)   =  NaN;
                SMRY.LatHalfMax_hit(snt,1)   =  NaN;
                SMRY.LatMaxPeak_hit(snt,1)   =  NaN;
                SMRY.LatFirstPeak_hit(snt,1) =  NaN;
                SMRY.LatPEP_hit(snt,1)       =  NaN;
                SMRY.LatTrise_hit(snt,1)     =  NaN;
            end
            
            try
                [~, latObj_right] = wz_spk_latency(spk.spiketimes(hit(right),:), latwin, sponttm, nBoot, smplsz,[],[],0);

                SMRY.LatFirstPeak_right(snt,1) =  latObj_right.first_peak;
                SMRY.LatMaxPeak_right(snt,1)   =  latObj_right.max_peak;
                SMRY.LatFirstResp_right(snt,1) =  latObj_right.peakresp(1);
                SMRY.LatMaxResp_right(snt,1)   =  latObj_right.peakresp(2);
                SMRY.LatHalfMax_right(snt,1)   =  latObj_right.halfmax;
                SMRY.LatMaxPeak_right(snt,1)   =  latObj_right.maxpeak;
                SMRY.LatFirstPeak_right(snt,1) =  latObj_right.frstpeak;
                SMRY.LatPEP_right(snt,1)       =  latObj_right.PEP;
                SMRY.LatTrise_right(snt,1)     =  latObj_right.Trise;
            catch
                SMRY.LatFirstPeak_right(snt,1) =  NaN;
                SMRY.LatMaxPeak_right(snt,1)   =  NaN;
                SMRY.LatFirstResp_right(snt,1) =  NaN;
                SMRY.LatMaxResp_right(snt,1)   =  NaN;
                SMRY.LatHalfMax_right(snt,1)   =  NaN;
                SMRY.LatMaxPeak_right(snt,1)   =  NaN;
                SMRY.LatFirstPeak_right(snt,1) =  NaN;
                SMRY.LatPEP_right(snt,1)       =  NaN;
                SMRY.LatTrise_right(snt,1)     =  NaN;
            end

            try
                [~, latObj_left] = wz_spk_latency(spk.spiketimes(hit(left),:), latwin, sponttm, nBoot, smplsz,[],[],0);

                SMRY.LatFirstPeak_left(snt,1) =  latObj_left.first_peak;
                SMRY.LatMaxPeak_left(snt,1)   =  latObj_left.max_peak;
                SMRY.LatFirstResp_left(snt,1) =  latObj_left.peakresp(1);
                SMRY.LatMaxResp_left(snt,1)   =  latObj_left.peakresp(2);
                SMRY.LatHalfMax_left(snt,1)   =  latObj_left.halfmax;
                SMRY.LatMaxPeak_left(snt,1)   =  latObj_left.maxpeak;
                SMRY.LatFirstPeak_left(snt,1) =  latObj_left.frstpeak;
                SMRY.LatPEP_left(snt,1)       =  latObj_left.PEP;
                SMRY.LatTrise_left(snt,1)     =  latObj_left.Trise;
            catch
                SMRY.LatFirstPeak_left(snt,1) =  NaN;
                SMRY.LatMaxPeak_left(snt,1)   =  NaN;
                SMRY.LatFirstResp_left(snt,1) =  NaN;
                SMRY.LatMaxResp_left(snt,1)   =  NaN;
                SMRY.LatHalfMax_left(snt,1)   =  NaN;
                SMRY.LatMaxPeak_left(snt,1)   =  NaN;
                SMRY.LatFirstPeak_left(snt,1) =  NaN;
                SMRY.LatPEP_left(snt,1)       =  NaN;
                SMRY.LatTrise_left(snt,1)     =  NaN;
            end
            
            try
                [~, latObj_pref] = wz_spk_latency(spk.spiketimes(hit(prefvis),:), latwin, sponttm, nBoot, smplsz,[],[],0);

                SMRY.LatFirstPeak_PrefVis(snt,1) =  latObj_pref.first_peak;
                SMRY.LatMaxPeak_PrefVis(snt,1)   =  latObj_pref.max_peak;
                SMRY.LatFirstResp_PrefVis(snt,1) =  latObj_pref.peakresp(1);
                SMRY.LatMaxResp_PrefVis(snt,1)   =  latObj_pref.peakresp(2);
                SMRY.LatHalfMax_PrefVis(snt,1)   =  latObj_pref.halfmax;
                SMRY.LatMaxPeak_PrefVis(snt,1)   =  latObj_pref.maxpeak;
                SMRY.LatFirstPeak_PrefVis(snt,1) =  latObj_pref.frstpeak;
                SMRY.LatPEP_PrefVis(snt,1)       =  latObj_pref.PEP;
                SMRY.LatTrise_PrefVis(snt,1)     =  latObj_pref.Trise;
            catch
                SMRY.LatFirstPeak_PrefVis(snt,1) =  NaN;
                SMRY.LatMaxPeak_PrefVis(snt,1)   =  NaN;
                SMRY.LatFirstResp_PrefVis(snt,1) =  NaN;
                SMRY.LatMaxResp_PrefVis(snt,1)   =  NaN;
                SMRY.LatHalfMax_PrefVis(snt,1)   =  NaN;
                SMRY.LatMaxPeak_PrefVis(snt,1)   =  NaN;
                SMRY.LatFirstPeak_PrefVis(snt,1) =  NaN;
                SMRY.LatPEP_PrefVis(snt,1)       =  NaN;
                SMRY.LatTrise_PrefVis(snt,1)     =  NaN;
            end

           %% target selection times
            resp_left  = nanmean(align_Clip.spikedensities(left,:));
            resp_right = nanmean(align_Clip.spikedensities(right,:));
            respdiff = abs(resp_left-resp_right);
            
            try
                Tdiff05  = wz_spk_RespDiff_ranksum(align_Clip.spikedensities(left,:), align_Clip.spikedensities(right,:), ...
                    [40, 200], 20, 0.05,  nBoot, smplsz, align_Clip.xtime);
                SMRY.TDT05(snt,1)  =  Tdiff05.DiffTime;
            catch
                SMRY.TDT05(snt,1)  =  NaN;
            end
            
            if(isfinite(SMRY.TDT05(snt,1)))               
                tdtpos = find(align_Clip.xtime == SMRY.TDT05(snt,1),1,'first');
                
                diffinit = respdiff(tdtpos);
                diff10   = respdiff(tdtpos+9);
                diff20   = respdiff(tdtpos+19);
                diff50   = respdiff(tdtpos+49);
                
                SMRY.TDT05_slope10(snt,1) = (diff10 - diffinit) / 10;
                SMRY.TDT05_area10(snt,1)  = sum(respdiff(tdtpos:tdtpos+9));
                
                SMRY.TDT05_slope20(snt,1) = (diff20 - diffinit) / 20;
                SMRY.TDT05_area20(snt,1)  = sum(respdiff(tdtpos:tdtpos+19));
                
                SMRY.TDT05_slope50(snt,1) = (diff50 - diffinit) / 50;                
                SMRY.TDT05_area50(snt,1)  = sum(respdiff(tdtpos:tdtpos+49));
            else
                SMRY.TDT05_slope10(snt,1) = NaN;
                SMRY.TDT05_area10(snt,1)  = NaN;
                SMRY.TDT05_slope20(snt,1) = NaN;
                SMRY.TDT05_area20(snt,1)  = NaN;
                SMRY.TDT05_slope50(snt,1) = NaN;                                
                SMRY.TDT05_area50(snt,1)  = NaN;                                
            end

            try
                Tdiff01  = wz_spk_RespDiff_ranksum(align_Clip.spikedensities(left,:), align_Clip.spikedensities(right,:), ...
                    [40, 200], 20, 0.01,  nBoot, smplsz, align_Clip.xtime);
                SMRY.TDT01(snt,1)  =  Tdiff01.DiffTime;
            catch
                SMRY.TDT01(snt,1)  =  NaN;
            end
            
            if(isfinite(SMRY.TDT01(snt,1)))               
                tdtpos = find(align_Clip.xtime == SMRY.TDT01(snt,1),1,'first');
                
                diffinit = respdiff(tdtpos);
                diff10   = respdiff(tdtpos+9);
                diff20   = respdiff(tdtpos+19);
                diff50   = respdiff(tdtpos+49);
                
                SMRY.TDT01_slope10(snt,1) = (diff10 - diffinit) / 10;
                SMRY.TDT01_area10(snt,1)  = sum(respdiff(tdtpos:tdtpos+9));
                
                SMRY.TDT01_slope20(snt,1) = (diff20 - diffinit) / 20;
                SMRY.TDT01_area20(snt,1)  = sum(respdiff(tdtpos:tdtpos+19));
                
                SMRY.TDT01_slope50(snt,1) = (diff50 - diffinit) / 50;                
                SMRY.TDT01_area50(snt,1)  = sum(respdiff(tdtpos:tdtpos+49));
            else
                SMRY.TDT01_slope10(snt,1) = NaN;
                SMRY.TDT01_area10(snt,1)  = NaN;
                SMRY.TDT01_slope20(snt,1) = NaN;
                SMRY.TDT01_area20(snt,1)  = NaN;
                SMRY.TDT01_slope50(snt,1) = NaN;                                
                SMRY.TDT01_area50(snt,1)  = NaN;                                
            end

            try
                Tdiff001 = wz_spk_RespDiff_ranksum(align_Clip.spikedensities(left,:), align_Clip.spikedensities(right,:), ...
                    [40, 200], 20, 0.001, nBoot, smplsz, align_Clip.xtime);
                SMRY.TDT001(snt,1) =  Tdiff001.DiffTime;
            catch
                SMRY.TDT001(snt,1)  =  NaN;
            end
            
            if(isfinite(SMRY.TDT001(snt,1)))               
                tdtpos = find(align_Clip.xtime == SMRY.TDT001(snt,1),1,'first');
                
                diffinit = respdiff(tdtpos);
                diff10   = respdiff(tdtpos+9);
                diff20   = respdiff(tdtpos+19);
                diff50   = respdiff(tdtpos+49);
                
                SMRY.TDT001_slope10(snt,1) = (diff10 - diffinit) / 10;
                SMRY.TDT001_area10(snt,1)  = sum(respdiff(tdtpos:tdtpos+9));
                
                SMRY.TDT001_slope20(snt,1) = (diff20 - diffinit) / 20;
                SMRY.TDT001_area20(snt,1)  = sum(respdiff(tdtpos:tdtpos+19));
                
                SMRY.TDT001_slope50(snt,1) = (diff50 - diffinit) / 50;                
                SMRY.TDT001_area50(snt,1)  = sum(respdiff(tdtpos:tdtpos+49));
            else
                SMRY.TDT001_slope10(snt,1) = NaN;
                SMRY.TDT001_area10(snt,1)  = NaN;
                SMRY.TDT001_slope20(snt,1) = NaN;
                SMRY.TDT001_area20(snt,1)  = NaN;
                SMRY.TDT001_slope50(snt,1) = NaN;                                
                SMRY.TDT001_area50(snt,1)  = NaN;                                
            end

            try
                Tdiff = wz_spk_RespDiff_diff(   align_Clip.spikedensities(left,:), align_Clip.spikedensities(right,:), ...
                    [40, 200], sponttm, 20, nBoot, smplsz, align_Clip.xtime);
                SMRY.Tdiff(snt,1)  =  Tdiff.DiffTime;
            catch
                SMRY.Tdiff(snt,1)  =  NaN;
            end
            
            if(isfinite(SMRY.Tdiff(snt,1)))               
                tdtpos = find(align_Clip.xtime == SMRY.Tdiff(snt,1),1,'first');
                
                diffinit = respdiff(tdtpos);
                diff10   = respdiff(tdtpos+9);
                diff20   = respdiff(tdtpos+19);
                diff50   = respdiff(tdtpos+49);
                
                SMRY.Tdiff_slope10(snt,1) = (diff10 - diffinit) / 10;
                SMRY.Tdiff_area10(snt,1)  = sum(respdiff(tdtpos:tdtpos+9));
                
                SMRY.Tdiff_slope20(snt,1) = (diff20 - diffinit) / 20;
                SMRY.Tdiff_area20(snt,1)  = sum(respdiff(tdtpos:tdtpos+19));
                
                SMRY.Tdiff_slope50(snt,1) = (diff50 - diffinit) / 50;                
                SMRY.Tdiff_area50(snt,1)  = sum(respdiff(tdtpos:tdtpos+49));
            else
                SMRY.Tdiff_slope10(snt,1) = NaN;
                SMRY.Tdiff_area10(snt,1)  = NaN;
                SMRY.Tdiff_slope20(snt,1) = NaN;
                SMRY.Tdiff_area20(snt,1)  = NaN;
                SMRY.Tdiff_slope50(snt,1) = NaN;                                
                SMRY.Tdiff_area50(snt,1)  = NaN;                                
            end

           %% ISI summary
            ISIobj = wz_spk_isi(spk.spiketimes(hit,:), sponttm, 0);

            SMRY.ISI_mean(snt,1)       = ISIobj.ISImean;
            SMRY.ISI_FF(snt,1)         = ISIobj.FF;
            SMRY.ISI_CV(snt,1)         = ISIobj.CV;
            SMRY.ISI_CV2(snt,1)        = ISIobj.CV2;
            SMRY.ISI_LV(snt,1)         = ISIobj.LV;
            SMRY.ISI_IR(snt,1)         = ISIobj.IR;
            SMRY.ISI_Bins2Peak(snt,1)  = ISIobj.RefRat2;
            SMRY.ISI_Bins2All(snt,1)   = ISIobj.RefRat1;
            SMRY.ISI_HalfMaxBin(snt,1) = ISIobj.HalfMaxBin;
            SMRY.logISI_mean(snt,1)    = ISIobj.logISI_mean;
            SMRY.logISI_median(snt,1)  = ISIobj.logISI_median;
            SMRY.logISI_CV(snt,1)      = ISIobj.logISI_CV;
            SMRY.logISI_skew(snt,1)    = ISIobj.logISI_skew;
            SMRY.logISI_kurtos(snt,1)  = ISIobj.logISI_kurtos;
            SMRY.logISI_hdt_dip(snt,1) = ISIobj.hdt_dip;
            SMRY.logISI_hdt_p(snt,1)   = ISIobj.hdt_p;
            SMRY.logISI_BC(snt,1)      = ISIobj.logISI_BC;

           %% plot summary
            if(maxrate < 2)
                maxrate = 10;
            end

           %% plot spike waveform
            subplot(plthdl(j,1));

            if(isfield(spk,'wave_file'))
                cwvfl = fullfile(rootdir, cSDir,'DSP',cDSp,spk.wave_file);
                if(exist(cwvfl,'file'))
                    cflwv = load(cwvfl);

                    [~,tpos] = intersect(cflwv.wave.spike_trials, hit);

                    if(~isempty(tpos))
                        SPKWV = wz_spk_get_wave(cflwv.wave.waves(tpos,:), cflwv.wave.thresh, 1, [], nfo, 1);
                        SMRY.WV_minmax(snt,1)   =  SPKWV.minmax;
                        SMRY.WV_trough(snt,1)   =  SPKWV.troughbin;
                        SMRY.WV_width(snt,1)    =  SPKWV.width;
                        SMRY.WV_duration(snt,1) =  SPKWV.duration;
                        SMRY.WV_area(snt,1)     =  SPKWV.area(1);
                        SMRY.WV_SNR(snt,1)      =  SPKWV.SNR;
                        SMRY.WV_SNR2(snt,1)     =  SPKWV.SNR2(1);

                        SMRY.WV_thr(snt,1)      = cflwv.wave.thresh;
                        SMRY.WV_SpkSign(snt,:)  = SPKWV.SpikeSign;
                    else
                        SMRY.WV_minmax(snt,1)   = NaN;
                        SMRY.WV_trough(snt,1)   = NaN;
                        SMRY.WV_width(snt,1)    = NaN;
                        SMRY.WV_duration(snt,1) = NaN;
                        SMRY.WV_area(snt,1)     = NaN;
                        SMRY.WV_SNR(snt,1)      = NaN;
                        SMRY.WV_SNR2(snt,1)     = NaN;

                        SMRY.WV_thr(snt,1)      = NaN;
                        SMRY.WV_SpkSign(snt,:)  = 'NA ';
                    end
                    % hline(cflwv.wave.thresh,'LineStyle',':');
                end
            end
            ylabel([char(sfx_arr(j)), detstr],'Interpreter','none', 'FontSize', 10);

            %% Plot SRT distributions
            subplot(plthdl(j,2));

            %             [h, sth]  = cdfplot(SRT_hit);
            %             set(h,'color','g', 'LineWidth', 2);
            hold on;

            if(~isempty(SRT_false))
                [h, stf]  = cdfplot(SRT_false);
                set(h,'color','k', 'LineWidth', 2);
                plot([stf.median, stf.median], [0, 0.5],'color','k', 'LineWidth', 2);
            end

            if(~isempty(SRT_left))
                [h, stl]  = cdfplot(SRT_left);
                set(h,'color','r', 'LineWidth', 2);
                plot([stl.median, stl.median], [0, 0.5],'color','r', 'LineWidth', 2);
            end

            if(~isempty(SRT_right))
                [h, str]  = cdfplot(SRT_right);
                set(h,'color','b', 'LineWidth', 2);
                plot([str.median, str.median], [0, 0.5],'color','b', 'LineWidth', 2);
            end

            xlim([min([SRT_hit(:); SRT_false(:)]), max([prctile(SRT_hit,90) ; prctile(SRT_false,90)])]);
            if(j==1)
                title('SRTs');
            else
                title([]);
            end

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);

            %% Stimulus onset
            subplot(plthdl(j,3));
            wz_spk_plot_hist(align_Clip, vis_plw, maxrate,0,10,0,100,1/3);
            
            if(isfinite(SMRY.LatHalfMax_hit(snt,1)))
                vline(SMRY.LatHalfMax_hit(snt,1), 'color','blue')
            end
            
            if(j==1)
                title('stimulus');
            end

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);

            %% Saccade onset
            subplot(plthdl(j,4));
            wz_spk_plot_hist(align_Sacc, sac_plw, maxrate,0,10,0,100,1/3);

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);
            if(j==1)
                title('saccade');
            end

            %% Reward
            subplot(plthdl(j,5));
            wz_spk_plot_hist(align_Rew, sac_plw, maxrate,0,10,0,100,1/3);

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);
            if(j==1)
                title('reward');
            end

            %% Error
            subplot(plthdl(j,6));
            if(~isempty(align_fFB.upperdense))
                try
                    wz_spk_plot_hist(align_fFB, sac_plw, maxrate,0,10,0,100,1/3);
                end
            end

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);
            if(j==1)
                title('error');
            end

            %% RespDiff Stim Onset
            subplot(plthdl(j,7));

            xvec = align_Clip.xtime;
            plpos = xvec >= vis_plw(1) & xvec <= vis_plw(2);

            vline(0,'color','k');
            vline(nanmedian(spk.Task.SRT(left)), 'color','r');
            vline(nanmedian(spk.Task.SRT(right)),'color','b');

            left_dense = align_Clip.spikedensities(left,:);
            Lmn  = nanmean(left_dense);
            Lste = nanstd(left_dense) ./ sqrt(sum(isfinite(left_dense)));
            shadedErrorBar(xvec,Lmn,Lste,'r',0)

            right_dense = align_Clip.spikedensities(right,:);
            Rmn  = nanmean(right_dense);
            Rste = nanstd(right_dense) ./ sqrt(sum(isfinite(right_dense)));
            shadedErrorBar(xvec,Rmn,Rste,'b',0)

            if(~isnan(SMRY.TDT05(snt,1)))
                vline(SMRY.TDT05(snt,1),'color','g','LineWidth',2);
            end

            maxY = max([max(Lmn(plpos)+Lste(plpos)),max(Rmn(plpos)+Rste(plpos))]);
            if(maxY < 4)
                maxY = 8;
            end
            ylim([0,maxY]);
            xlim(vis_plw);

            if(j==1)
                title('L vs R Stim');
            end

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);

            %% RespDiff Saccade Onset
            subplot(plthdl(j,8));

            xvec = align_Sacc.xtime;
            plpos = xvec >= sac_plw(1) & xvec <= sac_plw(2);

            vline(0,'color','k');

            left_dense = align_Sacc.spikedensities(left,:);
            Lmn  = nanmean(left_dense);
            Lste = nanstd(left_dense) ./ sqrt(sum(isfinite(left_dense)));
            shadedErrorBar(xvec,Lmn,Lste,'r',0)

            right_dense = align_Sacc.spikedensities(right,:);
            Rmn  = nanmean(right_dense);
            Rste = nanstd(right_dense) ./ sqrt(sum(isfinite(right_dense)));
            shadedErrorBar(xvec,Rmn,Rste,'b',0)

            maxY = max([max(Lmn(plpos)+Lste(plpos)),max(Rmn(plpos)+Rste(plpos))]);
            if(maxY < 4)
                maxY = 8;
            end
            ylim([0,maxY]);
            xlim(sac_plw);

            if(j==1)
                title('L vs R Sacc');
            end

            if(j==num_rows)
                xlabel('time [ms]');
            else
                xlabel([]);
            end
            ylabel([]);

            %% get oritation tuning

            %get direction statistics
            vis_peak_mu   = circ_mean(thv, peak_stim, 2);
            vis_peak_tune =    circ_r(thv, peak_stim./max(peak_stim), [], 2);

            vis_mean_mu   = circ_mean(thv, mean_stim./max(mean_stim),2);
            vis_mean_tune =    circ_r(thv, mean_stim./max(mean_stim), [], 2);

            mov_peak_mu   = circ_mean(thv, peak_sacc, 2);
            mov_peak_tune =    circ_r(thv, peak_sacc./max(peak_sacc), [], 2);

            mov_mean_mu   = circ_mean(thv, mean_sacc./max(mean_sacc),2);
            mov_mean_tune =    circ_r(thv, mean_sacc./max(mean_sacc), [], 2);

            SMRY.PrefOri_vis_peak(snt,1) = rad2deg(vis_peak_mu);
            SMRY.Tune_vis_peak(snt,1)    = vis_peak_tune;
            SMRY.PrefOri_vis_mean(snt,1) = rad2deg(vis_mean_mu);
            SMRY.Tune_vis_mean(snt,1)    = vis_mean_tune;
            SMRY.PrefOri_mov_peak(snt,1) = rad2deg(mov_peak_mu);
            SMRY.Tune_mov_peak(snt,1)    = mov_peak_tune;
            SMRY.PrefOri_mov_mean(snt,1) = rad2deg(mov_mean_mu);
            SMRY.Tune_mov_mean(snt,1)    = mov_mean_tune;

            % for stim onset
            subplot(plthdl(j,9));
            hold off
            wz_polar(trgtloc, peak_stim, '-r');
            hold on;
            h=polar(thv, peak_stim, 'or');
            set(h, 'MarkerFaceColor', 'r')

            wz_polar(trgtloc, mean_stim, '-b');
            h=polar(thv, mean_stim, 'ob');
            set(h, 'MarkerFaceColor', 'b')

            wz_polar(trgtloc, spontrrate, '-g');
            wz_polar(trgtloc, spontrrate+spontrerr, '--g');

            if(j==1)
                title('Vis Tune');
            end

            % get oritation tuning for saccade
            subplot(plthdl(j,10));
            hold off

            wz_polar(trgtloc, peak_sacc, '-r');
            hold on;
            h=polar(thv, peak_sacc, 'or');
            set(h, 'MarkerFaceColor', 'r')

            wz_polar(trgtloc, mean_sacc, '-b');
            h=polar(thv, mean_sacc, 'ob');
            set(h, 'MarkerFaceColor', 'b')

            wz_polar(trgtloc, spontrrate, '-g');
            wz_polar(trgtloc, spontrrate+spontrerr, '--g');

            if(j==1)
                title('Mov Tune');
            end
        end   % for(j=1:num_rows)

        if(do_plot == 1)
            hgexport(gcf, fullfile(odir,[cSess,'_',cDSp,'_overview.png']), hgexport('factorystyle'), 'Format', 'png');
            hgexport(gcf, fullfile(odir,[cSess,'_',cDSp,'_overview.eps']), hgexport('factorystyle'), 'Format', 'eps');
            close(fhndl);
        end
    end   %for(i=1:num_fls)
end   % for(sess=1:length(SessDir))

save(fullfile(odir, 'summary_table.mat'),'-struct', 'SMRY', '-v7');
TBL = struct2table(SMRY);
writetable(TBL,fullfile(odir, 'summary_table.dat'));

