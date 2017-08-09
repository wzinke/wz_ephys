function PLX2DB(fllst, srcpath, oroot, EV, keep_sess)
% Read in a plexon data file and save data it in a hierarchical structure.
%
% DESCRIPTION
%
%
%
% SYNTAX
%
%
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 30-Jan-2015 by wolf zinke
%
% ToDo: - Add option to update or overwrite current data base. Also, check
%         for for the file version that created the data, and in case it is
%         outdated re-create the data to ensure consistency.
%       - write summary/meta information into the session root directory
%         that identifies available data.



% ____________________________________________________________________________ %
%% get information about the current file version
% if(exist('wz_get_git_version','file'))
%     [gitversion, ~, GitInfo]    = wz_get_git_version(mfilename('fullpath'));
%     plxMat.TranslateFileVersion = gitversion;
%     plxMat.TranslateFileDate    = GitInfo.date;
% end

% ____________________________________________________________________________ %
%% check input data

% specify the file list
if(~exist('fllst','var') || isempty(fllst))
    [FileName,PathName] = uigetfile({'*.xls;*.xlsx;*.ods;*.csv;*.dat;*.txt;*.gnumeric'}, 'file list');
    fllst = fullfile(PathName,FileName);
end

% set source path where the plexon files are locates
if(~exist('srcpath','var') || isempty(srcpath))
    srcpath = uigetdir(cd, 'directory with plexon file');
end

% get the root path to the output. Note, that this is the parent directory of
% the individual animal.
if(~exist('oroot','var') || isempty(oroot))
    oroot = uigetdir(cd, 'root directory to safe converted files');
end

if(~exist(oroot,'dir'))
    mkdir(oroot);
end

if(~exist('keep_sess','var') || isempty(keep_sess))
    keep_sess = 0; % [Hz] - keep times at this time resolution, but only use integer numbers for space reasons
end

% ____________________________________________________________________________ %
%% get event code information
if(~exist('EV','var') || isempty(EV))
    [FileName,PathName] = uigetfile({'*.pro;*.mat;*.m'},'Event code definition file');
    EV = fullfile(PathName,FileName);
end

if(~isstruct(EV))
     [~,~,cflext] = fileparts(EV);   
    
     if(strcmpi(cflext,'.pro'))
        [EV] = TEMPO_get_eventdefs(EV);
     elseif(strcmpi(cflext,'.mat'))
         EV = load(EV);
     elseif(strcmpi(cflext,'.m'))
         [evpth, evstem] = fileparts(EV);
         if(~exist(evstem,'file'))  % stupid matlab can not run functions unless they are in the known paths...
             addpath(evpth);
         end
         eval(sprintf('EV = %s;', evstem)); 
     else
         error('Event file format not recognized!');
     end
end

% ____________________________________________________________________________ %
%% read in the file list

tbl = readtable(fllst);

[flpath, flstem] = fileparts(fllst);

%% check file

if(isempty(strmatch('Filename',char(fieldnames(tbl)))))
    error('File list requires a column with the name >Filename<!');
end

% create session list - account for the possibilitiy that there are more
% sessions on the same date.
sessIDs = [char(tbl.Date), char(tbl.session)];

% ____________________________________________________________________________ %
%% create required directories
NHPlst = cell2mat(unique(tbl.NHP));

for(i=1:length(NHPlst))
    switch NHPlst(i)
        case 'G'
            NHPname = 'Gauss';
        case 'H'
            NHPname = 'Helmholtz';
        case 'S'
            NHPname = 'Seymour';
        case 'Q'
            NHPname = 'Quincy';
        case 'D'
            NHPname = 'Darwin';
        case 'E'
            NHPname = 'Euler';
    end

    pmkdir(fullfile(oroot, NHPname));

    sessions = unique(sessIDs,'rows');

    for(p=1:size(sessions,1))
        cSess = sessions(p,:);

        pmkdir(fullfile(oroot, NHPname,cSess));

        pmkdir(fullfile(oroot, NHPname,cSess,'DSP'));
        pmkdir(fullfile(oroot, NHPname,cSess,'LFP'));
        pmkdir(fullfile(oroot, NHPname,cSess,'EEG'));
        pmkdir(fullfile(oroot, NHPname,cSess,'BHV'));
    end
end

% ____________________________________________________________________________ %
%% go through the file list
numfls=length(tbl.Filename);

dsp_cnt = 0;
lfp_cnt = 0;
eeg_cnt = 0;

for(i=1:numfls)

    NHP = char(tbl.NHP(i));
    rig = tbl.rig(i);

    switch NHP
        case 'G'
            NHPname = 'Gauss';
        case 'H'
            NHPname = 'Helmholtz';
    end

    if(strcmp(tbl.Specifier(i),'NA'))
        specstr = '';
    else
        specstr = ['_',char(tbl.Specifier(i))];
    end

    if(strcmp(tbl.fileSFX(i),'NA'))
        sfxstr = '';
    else
        sfxstr = ['_',char(tbl.fileSFX(i))];
    end

    cfl = char(tbl.Filename(i));

    if(~exist(fullfile(srcpath,cfl),'file'))
        warning(['File not found (',cfl,') -> skipped!']);
        continue;
    end

    disp(['processing ',cfl]);

    % get the plexon file in a pre-processed format
    plxMat = PLX_get_paradigm(fullfile(srcpath,cfl), EV, rig, NHP, keep_sess, 0);

    if(isempty(plxMat))
        warning('Not sufficient trials!');
        continue;
    end

    % loop over units
    disp(' ... processing DSP data:');
    for(u=1:plxMat.num_units)
        dsp_cnt = dsp_cnt + 1;
        dspname = char(plxMat.DSP_names(u));

        unitfile = [sessIDs(i,:),'_',dspname,'_',char(tbl.Paradigm(i)),sfxstr,specstr,'.mat'];
        disp(['     ... working on ',dspname]);

        DSP.NHP(dsp_cnt,1)          = tbl.NHP(i);
        DSP.Date(dsp_cnt,1)         = tbl.Date(i);
        DSP.Paradigm(dsp_cnt,1)     = tbl.Paradigm(i);
        DSP.Spec(dsp_cnt,1)         = tbl.Specifier(i);
        DSP.sessID(dsp_cnt,1)       = {sessIDs(i,:)};
        DSP.area(dsp_cnt,1)         = tbl.Area(i);
        DSP.File(dsp_cnt,1)         = {cfl};
        DSP.fileSFX(dsp_cnt,1)      = tbl.fileSFX(i);
        DSP.unitFile(dsp_cnt,1)     = {unitfile};
        DSP.Electrode(dsp_cnt,1)    = {plxMat.ElectrodeType};
        DSP.ap_coor(dsp_cnt,1)      = tbl.AP(i);
        DSP.ml_coor(dsp_cnt,1)      = tbl.ML(i);
        DSP.depth_uncorr(dsp_cnt,1) = tbl.Depth(i);
        DSP.DSPname(dsp_cnt,1)      = {dspname};
        DSP.ID(dsp_cnt,1)           = {[sessIDs(i,:),'_',dspname]};
        DSP.channel(dsp_cnt,1)      = str2num(dspname(4:5));
        DSP.neuron(dsp_cnt,1)       = char(dspname(6));
        DSP.NTrials(dsp_cnt,1)      = plxMat.NTrials;
        DSP.NCorrect(dsp_cnt,1)     = plxMat.NCorrect;

        unit_struct = [];
        unit_struct.NHP          = char(tbl.NHP(i));
        unit_struct.Date         = char(tbl.Date(i));
        unit_struct.Paradigm     = char(tbl.Paradigm(i));
        unit_struct.Spec         = char(tbl.Specifier(i));
        unit_struct.sessID       = sessIDs(i,:);
        unit_struct.area         = char(tbl.Area(i));
        unit_struct.File         = cfl;
        unit_struct.fileSFX      = char(tbl.fileSFX(i));
        unit_struct.Electrode    = char(plxMat.ElectrodeType);
        unit_struct.ap_coor      = char(tbl.AP(i));
        unit_struct.ml_coor      = char(tbl.ML(i));
        unit_struct.depth_uncorr = tbl.Depth(i);
        unit_struct.DSPname      = dspname;
        unit_struct.ID           = [sessIDs(i,:),'_',dspname];
        unit_struct.channel      = str2num(dspname(4:5));
        unit_struct.neuron       = char(dspname(6));

        unit_struct.Ntrials    = plxMat.NTrials;
        unit_struct.spiketimes = plxMat.DSP.(dspname);
        unit_struct.Task       = plxMat.Task;
        unit_struct.photodiode = plxMat.TrialMat_PDtimes;
        unit_struct.Infocodes  = plxMat.TrialMat_INFOcodes;
        unit_struct.EVcodes    = plxMat.TrialMat_EVcodes;
        unit_struct.EVtimes    = plxMat.TrialMat_EVtimes;
        unit_struct.SessTime   = plxMat.SessTime;

        uodir = fullfile(oroot, NHPname,sessIDs(i,:),'DSP',dspname);
        pmkdir(uodir);

        if(isfield(plxMat.DSP,[dspname,'_wave']))

            pmkdir(fullfile(uodir,'waves'));
            
            cwave = [];
            cwave.NHP       = unit_struct.NHP;
            cwave.Date      = unit_struct.Date;
            cwave.Paradigm  = unit_struct.Paradigm;
            cwave.Spec      = unit_struct.Spec;
            cwave.sessID    = unit_struct.sessID;
            cwave.area      = unit_struct.area;
            cwave.File      = unit_struct.File;
            cwave.fileSFX   = unit_struct.fileSFX;
            cwave.Electrode = unit_struct.Electrode;
            cwave.ap_coor   = unit_struct.ap_coor;
            cwave.ml_coor   = unit_struct.ml_coor;
            cwave.DSPname   = unit_struct.DSPname;
            cwave.ID        = unit_struct.ID;
            cwave.channel   = unit_struct.channel;
            cwave.neuron    = unit_struct.neuron;
            cwave.Ntrials   = unit_struct.Ntrials;
            cwave.wave      = plxMat.DSP.([dspname,'_wave']);

            [~, wvfile] = fileparts(unitfile);
            wvfile = fullfile('waves', [wvfile,'_waves.mat']);
            save(fullfile(uodir,wvfile), '-struct', 'cwave', '-v7');
            unit_struct.wave_file = wvfile;
        end

        unitfile = fullfile(uodir, unitfile);
        save(unitfile, '-struct', 'unit_struct', '-v7');
    end

    if(isfield(plxMat,'LFP'))
    % loop over LFP channels
        disp(' ... processing LFP data:');
        for(l=1:plxMat.num_lfp)            
            lfpname = char(plxMat.LFP_names(l));
            if(isfield(plxMat.LFP,lfpname))
                lfp_cnt = lfp_cnt + 1;
                chanfl = [sessIDs(i,:),'_',lfpname,'_',char(tbl.Paradigm(i)),sfxstr,specstr,'.mat'];
                disp(['     ... working on ',lfpname]);

                LFP.NHP(lfp_cnt,1)          = tbl.NHP(i);
                LFP.Date(lfp_cnt,1)         = tbl.Date(i);
                LFP.Paradigm(lfp_cnt,1)     = tbl.Paradigm(i);
                LFP.Spec(lfp_cnt,1)         = tbl.Specifier(i);
                LFP.sessID(lfp_cnt,1)       = {sessIDs(i,:)};
                LFP.area(lfp_cnt,1)         = tbl.Area(i);
                LFP.File(lfp_cnt,1)         = {cfl};
                LFP.fileSFX(lfp_cnt,1)      = tbl.fileSFX(i);
                LFP.chanFile(lfp_cnt,1)     = {chanfl};
                LFP.Electrode(lfp_cnt,1)    = {plxMat.ElectrodeType};
                LFP.ap_coor(lfp_cnt,1)      = tbl.AP(i);
                LFP.ml_coor(lfp_cnt,1)      = tbl.ML(i);
                LFP.depth_uncorr(lfp_cnt,1) = tbl.Depth(i);
                LFP.ID(lfp_cnt,1)           = {[sessIDs(i,:),'_', lfpname]};
                LFP.channel(lfp_cnt,1)      = plxMat.LFPchan(l);
                LFP.NTrials(lfp_cnt,1)      = plxMat.NTrials;
                LFP.NCorrect(lfp_cnt,1)     = plxMat.NCorrect;

                lfp_struct = [];
                lfp_struct.NHP          = char(tbl.NHP(i));
                lfp_struct.Date         = char(tbl.Date(i));
                lfp_struct.Paradigm     = char(tbl.Paradigm(i));
                lfp_struct.Spec         = char(tbl.Specifier(i));
                lfp_struct.sessID       = sessIDs(i,:);
                lfp_struct.area         = char(tbl.Area(i));
                lfp_struct.File         = cfl;
                lfp_struct.fileSFX      = char(tbl.fileSFX(i));
                lfp_struct.Electrode    = plxMat.ElectrodeType;
                lfp_struct.ap_coor      = char(tbl.AP(i));
                lfp_struct.ml_coor      = char(tbl.ML(i));
                lfp_struct.depth_uncorr = tbl.Depth(i);
                lfp_struct.ID           = [sessIDs(i,:),'_', lfpname];
                lfp_struct.channel      = plxMat.LFPchan(l);

                lfp_struct.Ntrials    = plxMat.NTrials;
                lfp_struct.timevec    = plxMat.timevec;
                lfp_struct.LFP        = plxMat.LFP.(lfpname);
        %         lfp_struct.LFP_z      = plxMat.LFP.([lfpname,'_z']);
                lfp_struct.info       = plxMat.LFP.([lfpname,'_info']);
                lfp_struct.Task       = plxMat.Task;
                lfp_struct.photodiode = plxMat.TrialMat_PDtimes;
                lfp_struct.Infocodes  = plxMat.TrialMat_INFOcodes;
                lfp_struct.EVcodes    = plxMat.TrialMat_EVcodes;
                lfp_struct.EVtimes    = plxMat.TrialMat_EVtimes;
                lfp_struct.SessTime   = plxMat.SessTime;

                lodir = fullfile(oroot, NHPname,sessIDs(i,:),'LFP',lfpname);
                pmkdir(lodir);
                lfpfile = fullfile(lodir,chanfl);
                save(lfpfile, '-struct', 'lfp_struct', '-v7');
            end
        end
        end

    if(isfield(plxMat,'EEG'))
    % loop over EEG channels
        disp(' ... processing EEG data:');
        for(e=1:length(plxMat.EEG_names))
            eeg_cnt = eeg_cnt + 1;
            eegname = char(plxMat.EEG_names(e));
            chanfl = [sessIDs(i,:),'_',eegname,'_',char(tbl.Paradigm(i)),sfxstr,specstr,'.mat'];
            disp(['     ... working on ',eegname]);

            EEG.NHP(eeg_cnt,1)          = tbl.NHP(i);
            EEG.Date(eeg_cnt,1)         = tbl.Date(i);
            EEG.Paradigm(eeg_cnt,1)     = tbl.Paradigm(i);
            EEG.Spec(eeg_cnt,1)         = tbl.Specifier(i);
            EEG.sessID(eeg_cnt,1)       = {sessIDs(i,:)};
            EEG.area(eeg_cnt,1)         = tbl.Area(i);
            EEG.File(eeg_cnt,1)         = {cfl};
            EEG.fileSFX(eeg_cnt,1)      = tbl.fileSFX(i);
            EEG.chanFile(eeg_cnt,1)     = {chanfl};
            EEG.Electrode(eeg_cnt,1)    = {plxMat.ElectrodeType};
            EEG.ap_coor(eeg_cnt,1)      = tbl.AP(i);
            EEG.ml_coor(eeg_cnt,1)      = tbl.ML(i);
            EEG.depth_uncorr(eeg_cnt,1) = tbl.Depth(i);
            EEG.ID(eeg_cnt,1)           = {[sessIDs(i,:),'_', eegname]};
            EEG.channel(eeg_cnt,1)      = plxMat.EEGchan(e);
            EEG.NTrials(eeg_cnt,1)      = plxMat.NTrials;
            EEG.NCorrect(eeg_cnt,1)     = plxMat.NCorrect;

            eeg_struct = [];
            eeg_struct.NHP          = char(tbl.NHP(i));
            eeg_struct.Date         = char(tbl.Paradigm(i));
            eeg_struct.Paradigm     = char(tbl.Paradigm(i));
            eeg_struct.Spec         = char(tbl.Specifier(i));
            eeg_struct.sessID       = sessIDs(i,:);
            eeg_struct.area         = char(tbl.Area(i));
            eeg_struct.File         = cfl;
            eeg_struct.fileSFX      = char(tbl.fileSFX(i));
            eeg_struct.Electrode    = plxMat.ElectrodeType;
            eeg_struct.ap_coor      = char(tbl.AP(i));
            eeg_struct.ml_coor      = char(tbl.ML(i));
            eeg_struct.depth_uncorr = tbl.Depth(i);
            eeg_struct.ID           = [sessIDs(i,:),'_', eegname];
            eeg_struct.channel      = plxMat.EEGchan(e);

            eeg_struct.Ntrials    = plxMat.NTrials;
            eeg_struct.timevec    = plxMat.timevec;
            eeg_struct.EEG        = plxMat.EEG.(eegname);
    %         eeg_struct.EEG_z      = plxMat.EEG.([eegname,'_z']);
            eeg_struct.info       = plxMat.EEG.([eegname,'_info']);
            eeg_struct.Task       = plxMat.Task;
            eeg_struct.photodiode = plxMat.TrialMat_PDtimes;
            eeg_struct.Infocodes  = plxMat.TrialMat_INFOcodes;
            eeg_struct.EVcodes    = plxMat.TrialMat_EVcodes;
            eeg_struct.EVtimes    = plxMat.TrialMat_EVtimes;
            eeg_struct.SessTime   = plxMat.SessTime;

            eodir = fullfile(oroot, NHPname,sessIDs(i,:),'EEG',eegname);
            pmkdir(eodir);
            eegfile = fullfile(eodir,chanfl);
            save(eegfile, '-struct', 'eeg_struct', '-v7');
        end
    end
    
    % save behaviour
    bhv_struct = [];
    bhv_struct.NHP     = char(tbl.NHP(i));
    bhv_struct.Date    = char(tbl.Paradigm(i));
    bhv_struct.sessID  = sessIDs(i,:);
    bhv_struct.area    = char(tbl.Area(i));
    bhv_struct.File    = cfl;
    bhv_struct.fileSFX = char(tbl.fileSFX(i));
    bhv_struct.Ntrials = plxMat.NTrials;
    bhv_struct.Task    = plxMat.Task;
    bhv_struct.timevec = plxMat.timevec;
    bhv_struct.Eye     = plxMat.Eye;

    bhvfile = fullfile(oroot, NHPname,sessIDs(i,:),'BHV',[sessIDs(i,:),'_',char(tbl.Paradigm(i)),sfxstr,specstr,'_behav.mat']);
    save(bhvfile, '-struct', 'bhv_struct', '-v7');

    % create a simple ascii table to be read into R as data frame
    bhv_tbl = [];
    bhv_tbl.NHP      = repmat(bhv_struct.NHP,bhv_struct.Ntrials,1);
    bhv_tbl.Date     = repmat(bhv_struct.Date,bhv_struct.Ntrials,1);
    bhv_tbl.sessID   = repmat(bhv_struct.sessID,bhv_struct.Ntrials,1);
    bhv_tbl.File     = repmat(bhv_struct.File,bhv_struct.Ntrials,1);
    bhv_tbl.fileSFX  = repmat(bhv_struct.fileSFX,bhv_struct.Ntrials,1);
    bhv_tbl.SessTime = plxMat.SessTime;

    tasknms = fieldnames(bhv_struct.Task);

    for(nm = 1:length(tasknms))
        cnm = char(tasknms(nm));
        if(length(bhv_struct.Task.(cnm)) == bhv_struct.Ntrials)
            bhv_tbl.(cnm) = bhv_struct.Task.(cnm)(:);
        end
    end

    BHV_TBL = struct2table(bhv_tbl);
    writetable(BHV_TBL, fullfile(oroot, NHPname,sessIDs(i,:),'BHV',[sessIDs(i,:),'_',char(tbl.Paradigm(i)),sfxstr,specstr,'_behav.dat']), 'Delimiter', ' ');
end

DSP_TBL = struct2table(DSP);
try
    writetable(DSP_TBL,fullfile(oroot, NHPname, [flstem,'_DSP.dat']), 'Delimiter', ' ');
    save(              fullfile(oroot, NHPname, [flstem,'_DSP.mat']),'-struct', 'DSP', '-v7');
end
try
    writetable(DSP_TBL,fullfile(flpath, [flstem,'_DSP.dat']), 'Delimiter', ' ');
    save(              fullfile(flpath, [flstem,'_DSP.mat']),'-struct', 'DSP', '-v7');
end

LFP_TBL = struct2table(LFP);
try
    writetable(LFP_TBL,fullfile(oroot, NHPname, [flstem,'_LFP.dat']), 'Delimiter', ' ');
    save(              fullfile(oroot, NHPname, [flstem,'_LFP.mat']),'-struct', 'LFP', '-v7');
end
try
    writetable(LFP_TBL,fullfile(flpath, [flstem,'_LFP.dat']), 'Delimiter', ' ');
    save(              fullfile(flpath, [flstem,'_LFP.mat']),'-struct', 'LFP', '-v7');
end

EEG_TBL = struct2table(EEG);
try
    writetable(EEG_TBL,fullfile(oroot, NHPname, [flstem,'_EEG.dat']), 'Delimiter', ' ');
    save(              fullfile(oroot, NHPname, [flstem,'_EEG.mat']),'-struct', 'EEG', '-v7');
end
try
    writetable(EEG_TBL,fullfile(flpath, [flstem,'_EEG.dat']), 'Delimiter', ' ');
    save(              fullfile(flpath, [flstem,'_EEG.mat']),'-struct', 'EEG', '-v7');
end
