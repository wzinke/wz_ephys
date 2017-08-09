function  DSP_plot_DSPmatrix(SessDir, odir)
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
%       SessDir         path to session directory which has to contain an 'LFP' subdirectory
%
%       odir            Output directory to save the summary plots (default is <SessDir>/DSP/plots)
%
%
%
%   Output:
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 15-Feb-2015 by wolf zinke
%
% ToDo:  - add loop to run over several input files
%        - write summary table with basic neuron descriptors

mincorrtrials = 40; % Number of trials for a condition to be plotted
maxspkrate    = 60; % maxinum firing rate during the pre-stim period to be accepted as usefull unit.

ktype = 'exp';
kwdth = 20;

vistune = [  75 125];
movtune = [ -50   0];
sponttm = [-250   0];

vis_win = [-500 1000];
sac_win = [-750 750];

% if not specified use GUI to get the file
if(~exist('SessDir','var') || isempty(SessDir))
    SessDir = uigetdir('path to session data');
end

if(~exist('odir','var') || isempty(odir))
    do_plot = 0;
else
    do_plot = 1;
    pmkdir(odir);
end

% ____________________________________________________________________________ %
%% get extension list
[~, cSess] = fileparts(SessDir);

DSPlst  = dir(fullfile(SessDir,'DSP','DSP*')); % get available DSP files
num_fls = length(DSPlst);

for(d=1:num_fls)
       
    cDSp = DSPlst(d).name;
    
    disp([char(cSess),'-',char(cDSp)]);
    
    % get available mat files
    mlst = dir(fullfile(SessDir,'DSP',cDSp,[cSess,'_',cDSp,'_*.mat']));
    numMat = length(mlst);
    
    % get a list of available paradigms
    dsp_arr   = [];
    task_arr  = {};
    setsz_arr = [];
    sfx_arr   = {};
    
    for(j=1:numMat)
        [~,cdata] = fileparts(mlst(j).name);
        splpos    = find(cdata == '_');
        fileSFX   = cdata(splpos(2)+1:end);
        
        cfl = fullfile(SessDir,'DSP',cDSp,[cSess,'_',cDSp,'_',fileSFX,'.mat']);
        
        dsp = load(cfl);
        
        % make a quick check if the firing rate is excessively high and discard if so
        spospk = sum(dsp.spiketimes(:) >= sponttm(1) & dsp.spiketimes(:) <= sponttm(2));
        spospkrt = 1000 * spospk / dsp.Ntrials /(diff(sponttm)+1);
        
        if(spospkrt > maxspkrate)
            continue;
        end
        
        p = dsp.Task.Correct == 1;
        
        tsks = unique(dsp.Task.TaskType);
        
        for(t=1:length(tsks))
            tp = strcmp(dsp.Task.TaskType, tsks(t)) & p == 1;
            
            if(sum(tp) > mincorrtrials)
                if(strcmp(tsks(t), 'MG'))
                    dsp_arr   = [dsp_arr; dsp];
                    task_arr  = [task_arr; 'MG'];
                    setsz_arr = [setsz_arr; NaN];
                    sfx_arr   = [sfx_arr; fileSFX];
                    
                elseif(strcmp(tsks(t), 'Cap'))
                    
                    if(sum(dsp.Task.Singleton(tp) == 1) > mincorrtrials)
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'Cap'];
                        setsz_arr = [setsz_arr; 0];
                        sfx_arr   = [sfx_arr; fileSFX];
                    end
                    
                    if(sum(dsp.Task.Singleton(tp) == 0) > mincorrtrials)
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'Cap'];
                        setsz_arr = [setsz_arr; 1];
                        sfx_arr   = [sfx_arr; fileSFX];
                    end                   
                    
                elseif(strcmp(tsks(t), 'Det'))
                    dsp_arr   = [dsp_arr; dsp];
                    task_arr  = [task_arr; 'Det'];
                    setsz_arr = [setsz_arr; 1];
                    sfx_arr   = [sfx_arr; fileSFX];
                    
                elseif(strcmp(tsks(t), 'Search'))
                     if(sum(dsp.Task.SetSize(tp) == 2) > mincorrtrials)
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'Search'];
                        setsz_arr = [setsz_arr; 2];
                        sfx_arr   = [sfx_arr; fileSFX];
                    end                   
                   
                     if(sum(dsp.Task.SetSize(tp) == 4) > mincorrtrials)
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'Search'];
                        setsz_arr = [setsz_arr; 4];
                        sfx_arr   = [sfx_arr; fileSFX];
                     end      
                    
                     if(sum(dsp.Task.SetSize(tp) == 8) > mincorrtrials)
                        dsp_arr   = [dsp_arr; dsp];
                        task_arr  = [task_arr; 'Search'];
                        setsz_arr = [setsz_arr; 8];
                        sfx_arr   = [sfx_arr; fileSFX];
                    end                   
                end
            end
        end        
    end
    
    if(isempty(dsp_arr))
        continue;
    end
    
    % ____________________________________________________________________________ %
    %% create plot outline:  Wave | Stim on | Sacc on | Reward | Error | LvsR Stim | LvsR Sacc | Ori Stim | Ori Sacc
    num_rows = length(dsp_arr);
    num_cols  = 9;   % change to 9 and add tuning plots for vis and mov
    
    figsz    = [0 0 num_cols*400 num_rows*260];
       
    fhndl  = figure('Name', [SessDir, '  -  ', cDSp], 'Position', figsz, 'Renderer', 'Painters');
    
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
  
    for(j=1:num_rows)
        disp(['...', char(sfx_arr(j)), ' - ', char(task_arr(j))]);

        spk = dsp_arr(j);
        
        % get relevant trial positions
        switch char(task_arr(j))
            case 'MG'
                tskp = strcmp(spk.Task.TaskType, task_arr(j)) & spk.Task.SRT > 120;
                detstr = [];
            case {'Search', 'Det'}
                tskp = strcmp(spk.Task.TaskType, task_arr(j)) & spk.Task.SRT > 120 & ...
                       spk.Task.SetSize ==  setsz_arr(j) & spk.Task.IsCatch == 0;
                if(setsz_arr(j) == 1)
                    detstr = [' - Detection'];
                else
                    detstr = [' - Set Size ',num2str(setsz_arr(j))];
                end
            case 'Cap'
                tskp = strcmp(spk.Task.TaskType, task_arr(j)) & spk.Task.SRT > 120 & ...
                       spk.Task.Singleton == setsz_arr(j) & spk.Task.IsCatch == 0;
                if(setsz_arr(j) == 1)
                    detstr = [' - With Singleton'];
                else
                    detstr = [' - No Singleton'];
                end
        end
        
        hit   = tskp == 1 & strcmp(spk.Task.TaskType, task_arr(j)) & spk.Task.Correct == 1; % select correct trials
        false = tskp == 1 & strcmp(spk.Task.TaskType, task_arr(j)) & spk.Task.error   == 1; % select false response
        left  = spk.Task.TargetLoc(hit) == 180;
        right = spk.Task.TargetLoc(hit) ==   0;
        
        align_Clip = wz_spk_density(wz_get_SPKobj(spk.spiketimes(hit,:),   vis_win, [], spk.Task.SRT(hit)),    ktype, kwdth);
        align_Sacc = wz_spk_density(wz_get_SPKobj(spk.spiketimes(hit,:),   sac_win,     spk.Task.SRT(hit)),    ktype, kwdth);
        align_Rew  = wz_spk_density(wz_get_SPKobj(spk.spiketimes(hit,:),   sac_win,     spk.Task.Reward(hit)), ktype, kwdth);
        align_fFB  = wz_spk_density(wz_get_SPKobj(spk.spiketimes(false,:), sac_win,     spk.Task.ErrorTone(false)), ktype, kwdth);

        vis_plw = [-100 ; prctile(spk.Task.SRT(hit),90)];
        sac_plw = [-1 *   prctile(spk.Task.SRT(hit),90); 100];

        stimpl = align_Clip.xtime > vis_plw(1) & align_Clip.xtime < vis_plw(2);
        sacmpl = align_fFB.xtime  > sac_plw(1) & align_fFB.xtime  < sac_plw(2);
                       
        if(~isempty(align_fFB.upperdense))
            maxrate = 1.1*max( [namax(align_Clip.upperdense(stimpl)); ...
                                namax(align_Sacc.upperdense(sacmpl)); ...
                                namax( align_Rew.upperdense(sacmpl)); ...
                                namax( align_fFB.upperdense(sacmpl))]);
        else
            maxrate = 1.1*max( [namax(align_Clip.upperdense(stimpl)); ...
                                namax(align_Sacc.upperdense(sacmpl)); ...
                                namax( align_Rew.upperdense(sacmpl))]);
        end
        
        if(maxrate < 2)
            maxrate = 10;
        end
        
        if(j==1)
            nfo = [cSess, '  -  ', cDSp];
        else
            nfo = [];
        end
        
       %% plot spike waveform
        subplot(plthdl(j,1));
        
        % quick and dirty fix (use intersect for a better solution!):
        % SCREW IT FOR NOW!
%         ctrials = find(hit==1);
%         
%         wvsel = [];
%         for(wt = 1 : length(ctrials))
%             wvsel = [wvsel; find(spk.wave.spike_trials == ctrials(wt))];
%         end
%         if(~isempty(wvsel))
        if(~isempty(spk.wave.waves))
                wz_spk_get_wave(spk.wave.waves, 1,[],nfo,1);
        end
        hline(spk.wave.thresh,'LineStyle',':');
%         end
        ylabel([char(sfx_arr(j)), detstr],'Interpreter','none', 'FontSize', 10);
        
       %% Stimulus onset
        subplot(plthdl(j,2));
        wz_spk_plot_hist(align_Clip, vis_plw, maxrate,0,10,0,100,1/3);
        if(j==1)
            title('stimulus');
        end
       
        xlabel([]); ylabel([]);

       %% Saccade onset
        subplot(plthdl(j,3));
        wz_spk_plot_hist(align_Sacc, sac_plw, maxrate,0,10,0,100,1/3);
       
        xlabel([]); ylabel([]);
        if(j==1)
            title('saccade');
        end
        
       %% Reward 
        subplot(plthdl(j,4));
        wz_spk_plot_hist(align_Rew, sac_plw, maxrate,0,10,0,100,1/3);

        xlabel([]); ylabel([]);
        if(j==1)
            title('reward');
        end
        
       %% Error
        subplot(plthdl(j,5));
        if(~isempty(align_fFB.upperdense))
            try
                wz_spk_plot_hist(align_fFB, sac_plw, maxrate,0,10,0,100,1/3);
            end
        end

        xlabel([]); ylabel([]);
        if(j==1)
            title('error');
        end
        
       %% RespDiff Stim Onset
        subplot(plthdl(j,6));
        
        xvec = align_Clip.xtime;
        plpos = xvec >= vis_plw(1) & xvec <= vis_plw(2);
        
        vline(0,'color','g');
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
        
        maxY = max([max(Lmn(plpos)+Lste(plpos)),max(Rmn(plpos)+Rste(plpos))]);
        if(maxY < 4)
            maxY = 8;
        end
        ylim([0,maxY]);
        xlim(vis_plw);

        if(j==1)
            title('L vs R Stim');
        end
        
       %% RespDiff Saccade Onset
        subplot(plthdl(j,7));

        xvec = align_Sacc.xtime;
        plpos = xvec >= sac_plw(1) & xvec <= sac_plw(2);

        vline(0,'color','g');
        
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
        
        %% get oritation tuning
        
        spntwin = align_Clip.xtime > sponttm(1) & align_Clip.xtime < sponttm(2);
        spontrrate = mean(nanmean(align_Clip.spikedensities(:,spntwin)));
        spontrerr  =  std(nanmean(align_Clip.spikedensities(:,spntwin),2));
        
        trgtloc = sort(unique(spk.Task.TargetLoc(hit)));
        
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
        
        if(j==1)
            title('Vis Tune');
        end
        
        % for stim onset
        subplot(plthdl(j,8));
        hold off
        wz_polar(trgtloc, peak_stim, '-r');
        hold on;
        h=polar(deg2rad(trgtloc), peak_stim', 'or');
        set(h, 'MarkerFaceColor', 'r')
        
        wz_polar(trgtloc, mean_stim, '-b');
        h=polar(deg2rad(trgtloc), mean_stim', 'ob');
        set(h, 'MarkerFaceColor', 'b')
        
        wz_polar(trgtloc, spontrrate, '-g');
        wz_polar(trgtloc, spontrrate+spontrerr, '--g');
        
        % get oritation tuning for saccade
        subplot(plthdl(j,9));
         hold off
               
        wz_polar(trgtloc, peak_sacc, '-r');
        hold on;
        h=polar(deg2rad(trgtloc), peak_sacc', 'or');
        set(h, 'MarkerFaceColor', 'r')
        
        wz_polar(trgtloc, mean_sacc, '-b');
        h=polar(deg2rad(trgtloc), mean_sacc', 'ob');
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

