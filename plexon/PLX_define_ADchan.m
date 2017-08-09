function [ADchan] = PLX_define_ADchan(rig, NHP, dtnum)
% Define the default layout of analog channels for a given rig.
%
% DESCRIPTION
%       This functions contains rig and animal specific assignments of analog
%       channels with some meaningful labels. The output is used to 'interprete'
%       and rename analog channels in the struct for the plexon data that was
%       created with PLX2Mat.
%
% SYNTAX
%
%   [plxdt] = PLX_define_ADchan(rig)
%
%   Input:
%
%       rig  - string that identifies the rig set up
%       NHP  - in case there are NHP specific adjustment in the configuration.
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 22-Jan-2015 by wolf zinke

switch rig
    case 28

        ADchan.Stim1  =  61;
        ADchan.Stim2  =  62;
        
        ADchan.EyeX  =  63;
        ADchan.EyeY  =  64;
        ADchan.Pupil =  [];

        switch upper(NHP)
            case 'D'
                ADchan.EEG.Fz  = 33;
                ADchan.EEG.F3  = 34;
                ADchan.EEG.F4  = 35;
                ADchan.EEG.Cz  = 36;
                ADchan.EEG.C3  = 37;
                ADchan.EEG.C4  = 38;
                ADchan.EEG.T5  = 39;
                ADchan.EEG.T6  = 40;
 
            case 'G'
                ADchan.EEG.Oz  = 33;
                ADchan.EEG.P4  = 34;
                ADchan.EEG.P3  = 35;
                ADchan.EEG.O2  = 36;
                ADchan.EEG.O1  = 37;
                ADchan.EEG.POz = 38;
                ADchan.EEG.Oz2 = 41;
                ADchan.EEG.OL  = 42;
                ADchan.EEG.OR  = 43;
                ADchan.EEG.PO3 = 44;
                ADchan.EEG.PO4 = 45;
                ADchan.EEG.Oz1 = 46;
                
                ADchan.EEG.US  = 57;
                ADchan.EEG.Reference = 40;

            case 'H'
                ADchan.EEG.FCz = 33;
                ADchan.EEG.P4  = 34;
                ADchan.EEG.P3  = 35;
                ADchan.EEG.O2  = 36;
                ADchan.EEG.O1  = 37;
                ADchan.EEG.POz = 38;
                ADchan.EEG.OR  = 41;
                ADchan.EEG.PO4 = 42;
                ADchan.EEG.Oz  = 43;
                ADchan.EEG.PO3 = 44;
                ADchan.EEG.OL  = 45;
                ADchan.EEG.C3  = 46;
                
                ADchan.EEG.US  = 57;
                ADchan.EEG.Reference = 40;

            otherwise
                error(['undefined data for ',NHP,' in rig ', rig]);
        end
        
    case 30
        if(dtnum >= datenum('1-Mar-2010'))
            ADchan.EyeX  =  47;
            ADchan.EyeY  =  48;
            ADchan.Pupil =  45;
            ADchan.Stim1  = 44;
        else
            ADchan.EyeX  =  18;
            ADchan.EyeY  =  19;
            ADchan.Pupil =  [];
        end
        
        warning('channel assignment needs to be implemented!');
        
    otherwise
        error(['No configuration specified for rig ', rig]);
end

% ____________________________________________________________________________ %
%% keep this the last step, not definition afterwards (except you know what you are doing)
def_fields = fieldnames(ADchan);

ADchan.NHP = NHP;
ADchan.rig = rig;

ADchan.is_eeg = [];
ADchan.ChanNum  = [];
ADchan.ChanName = {};

for (f = 1:length(def_fields))

    cfnm = char(def_fields(f));

    if(strncmp(cfnm,'EEG',3))
        more_def_fields = fieldnames(ADchan.EEG);
        for (e = 1:length(more_def_fields))
            cfnm = char(more_def_fields(e));
            cval = ADchan.EEG.(cfnm);
            if(~isempty(cval))
                ADchan.ChanNum  = [ADchan.ChanNum; cval];
                ADchan.ChanName = [ADchan.ChanName; cfnm];
                ADchan.is_eeg = [ADchan.is_eeg; 1];
            end
        end
    else
        cval = ADchan.(cfnm);
        if(~isempty(cval))
            ADchan.ChanNum  = [ADchan.ChanNum; cval];
            ADchan.ChanName = [ADchan.ChanName; cfnm];
            ADchan.is_eeg = [ADchan.is_eeg; 0];
        end
    end
end





