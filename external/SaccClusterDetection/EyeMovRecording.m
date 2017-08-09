classdef EyeMovRecording < ClusterDetection.DataDB
    %EYEMOVEMENTRECORDING Parent class for all recordings containing eye
    %movements
    %   Objects of this class represent single eye movement recordings
    %
    % Jorge Otero-Millan, jorgeoteromillan@gmail.com 2/17/2014
    %
    
    
    properties (Constant = true)
        
        % labels for the columns of position, velocity and acceleration
        % matrix
        
        LH = 1; % Left eye horizontal
        LV = 2; % Left eye vertical
        RH = 3; % Right eye horizontal
        RV = 4; % Right eye vertical
        LP = 5; % Left eye polar
        RP = 6; % Right eye polar
        
        % labels for the columns of the trial matrix
        TRIALSTART = 1;
        TRIALSTOP = 2;
        
        TRIAL_CROP_BEGINING = 1; % seconds to ignore at the begining of the trial
        TRIAL_CROP_END = 0.01; % seconds to ignore at the end of the trial
        TRIAL_MINIMUM_DURATION = 1; % seconds
        
        VALID_MINIMUM_DURATION = 0.1; % discard periods of continuos good data that are too short
        
    end
    
    properties % stored in the object (memory) SMALL VARIABLES
        samplerate % samples per second
        nsamples % total number of samples
        
        hasLeftEye
        hasRightEye
        trials
    end
    
    properties( Dependent = true ) % not stored in the object (memory) BIG VARIABLES
        
        time % timestamps of the recording
        
        position % position data
        velocity % velocity data
        acceleration % acceleration data
        
        valid % flags for each sample, 1 valid recording, 0 not valid (out of trial or blinks)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get functions for dependent variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function time = get.time(this)
            samples = this.ReadVariable('samples');
            time = samples(:,1);
        end
        
        function position = get.position(this)
            samples = this.ReadVariable('samples');
            position = samples(:,2:5);
        end
        
        function posAvg = getAvgPositionAcrossEyes( this )
            pos = this.position;
            if ( this.hasLeftEye && this.hasRightEye )
                posAvg = [mean(pos(:,[this.LH this.RH]),2) mean(pos(:,[this.LV this.RV]),2)];
            elseif ( this.hasLeftEye )
                posAvg = pos(:,[this.LH this.LV]);
            elseif ( this.hasRightEye )
                posAvg = pos(:,[this.RH this.RV]);
            end
        end
        
        function velocity = get.velocity(this)
            pos = this.position;
            velocity = zeros(length(pos(:,1)),6);
            velocity(:,[this.LH this.LV this.RH this.RV]) = this.Differenciate( pos(:,[this.LH this.LV this.RH this.RV]), this.samplerate );
            velocity(:,[this.LP this.RP]) = this.CartToPolar( velocity(:,[this.LH this.RH]), velocity(:,[this.LV this.RV]) );
        end
        
        function v = getAvgPolarVelocityAcrossEyes( this )
            pos = this.position;
            if ( this.hasLeftEye && this.hasRightEye )
                posh = mean([pos(:,this.LH ),pos(:,this.RH)],2);
                posv = mean([pos(:,this.LV ),pos(:,this.RV)],2);
                pos = [posh posv];
            elseif ( this.hasLeftEye )
                pos = pos(:,[this.LH this.LV]);
            elseif ( this.hasRightEye )
                pos = pos(:,[this.RH this.RV]);
            end
            
            vel = this.Differenciate( pos, this.samplerate );
            v = this.CartToPolar( vel(:,1), vel(:,2) );
        end
        
        function acceleration = get.acceleration(this)
            vel = this.velocity;
            acceleration = zeros(length(vel(:,1)),6);
            acceleration(:,[this.LH this.LV this.RH this.RV]) = this.Differenciate( vel(:,[this.LH this.LV this.RH this.RV]), this.samplerate );
            acceleration(:,[this.LP this.RP]) = this.CartToPolar( acceleration(:,[this.LH this.RH]), acceleration(:,[this.LV this.RV]) );
        end
        
        function a = getAvgPolarAccelerationAcrossEyes( this )
            acceleration = this.acceleration;
            if ( this.hasLeftEye && this.hasRightEye )
                a = mean([acceleration(:,this.LP),acceleration(:,this.RP)],2);
            elseif ( this.hasLeftEye )
                a = acceleration(:,this.LP);
            elseif ( this.hasRightEye )
                a = acceleration(:,this.RP);
            end
        end
        
        function valid = get.valid( this )
            blinkYesNo = this.ReadVariable('blinkYesNo');
            trialYesNo = this.CreateYesNo( this.trials(:,this.TRIALSTART), this.trials(:,this.TRIALSTOP), size(blinkYesNo,1) );
            
            % valid samples are samples within a trial and not in a blink
            valid = ~blinkYesNo & trialYesNo;
            
            % remove chunks of valid data that are too short 
            validStart = find(diff([0;valid])>0);
            validEnd = find(diff([valid;0])<0);
            validDuration = (validEnd-validStart)./this.samplerate;
            goodChunks = find(validDuration >  this.VALID_MINIMUM_DURATION);
            validChunks = [validStart(goodChunks) validEnd(goodChunks)];
            
            valid = this.CreateYesNo( validChunks(:,1), validChunks(:,2), size(blinkYesNo,1) );
        end
    end
    
    methods
        
        function this = Init( this )
            info = this.ReadVariable('info');
            this.nsamples = info.nSamples;
            this.samplerate	= info.samplerate;
            
            this.trials = info.trials;
            
            this.hasLeftEye = info.hasLeftEye;
            this.hasRightEye = info.hasRightEye;
        end
        
        function trials = GetTrials( this )
            
            time = this.time;
            % trials
            isInTrial = [0;diff(time)]==round(median(diff(time))); % jumps in timestamps
            
            trialStartIdx = find(diff([0;isInTrial])>0);
            trialEndIdx = find(diff([isInTrial;0])<0);
                        
            % adjust the begining and end of trials
            trialStartIdx = trialStartIdx + this.TRIAL_CROP_BEGINING*this.samplerate;
            trialEndIdx = trialEndIdx - this.TRIAL_CROP_END*this.samplerate;
            
            % select long enough trials
            trialDuration = (trialEndIdx-trialStartIdx)./this.samplerate;
            goodTrials = find(trialDuration >  this.TRIAL_MINIMUM_DURATION);
            trials = [trialStartIdx(goodTrials) trialEndIdx(goodTrials)];
            
        end
        
        function [saccades stats] = FindSaccades( this)
            
            saccadeDetector = ClusterDetection.SaccadeDetector.getSaccadeDetector( );
            [saccades stats ] = saccadeDetector.FindSaccades( this );   
            
            this.WriteVariable( saccades, 'saccades');
            this.WriteVariable( stats, 'stats');
        end
        
        function v = Differenciate( this, x, samplerate )
            N = length(x(:,1));            % length of the time series
            M = length(x(1,:));
            
            w = window(@bartlett,6);
            xf = filter(w/sum(w),1,x);
            dxf = diff(xf)*samplerate;
            
            v = zeros(N,M);
            v(1:end-3,:) = dxf(3:end,:);
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Static methods to create new recording objects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods( Static = true )
        
        % Experiment object factory from already created recording in DB
        function recording = Load( folder, session )
            
            % create Recording object
            recording = ClusterDetection.EyeMovRecording();
            
            % init DB and object
            recording.InitDB( folder, session );
            recording.Init();
        end
        
        % Experiment object factory from simple variables
        function recording = Create( folder, session, samples, blinkYesNo, samplerate, importInfo )
            try
                
                recording = ClusterDetection.EyeMovRecording();
                recording.InitDB( folder, session );
                
                if ( exist( 'importInfo','var') )
                    info.import = importInfo;
                end
                
                info.nSamples = size(samples,1);
                info.samplerate = samplerate;
                recording.samplerate = info.samplerate;
                
                info.hasLeftEye = (sum(~isnan( samples(:,2) )) > 0); %TODO more ways to determine if eyes are present ...
                info.hasRightEye = (sum(~isnan( samples(:,4) )) > 0);
                
                recording.WriteVariable( samples, 'samples');
                
                % save data in DB in new session
                recording.WriteVariable( blinkYesNo, 'blinkYesNo');
                
                info.trials = recording.GetTrials( );
                recording.WriteVariable( info, 'info');
                
                recording.Init();
                
            catch me
                rethrow(me)
            end
        end
        
    end
    
    
    
    %% Utility functions
    methods ( Static = true );
        
        function p = CartToPolar( x, y )
            p = sqrt( x.^2 + y.^2 );
        end
        
        function yesNo = CreateYesNo( starts, stops, length )
            % TODO: check for errors in starts and stops
            
            yesNo = zeros(length,1);
            yesNo(starts) = 1;
            yesNo(stops) = -1;
            yesNo = cumsum(yesNo);
        end
        
        
        function [lb rb lm rm] = FindBinocularEvents( l, r )
            ls = l(:,1);
            le = l(:,2);
            rs = r(:,1);
            re = r(:,2);
            
            binocs2 = [];
            
            % for each microsaccade in the left eye
            for i=1:length(ls)
                e = ones(size(re))*le(i);
                s = ones(size(rs))*ls(i);
                row = max(min( re, e) - max(rs, s),0);
                % row has the overlap of this left microsaccade with all the right
                % microsaccades
                if( sum(row) > 0 )
                    % we get the right microsaccade that overlaps more with the
                    % left one and we check that the left one is also the one that
                    % overlaps more with the right one.
                    [v, max_index_r] = max( row );
                    e = ones(size(le))*re(max_index_r);
                    s = ones(size(ls))*rs(max_index_r);
                    col = max( min( le, e) - max(ls, s), 0 );
                    [v, max_index_l] = max( col );
                    if (max_index_l == i )
                        binocs2(end+1, : ) = [max_index_l;max_index_r];
                    end
                end
            end
            if ( ~isempty( binocs2 ) )
                lb = binocs2(:,1);
                rb = binocs2(:,2);
            else
                lb = [];
                rb = [];
            end
            lm = setdiff(1:length(l),lb);
            rm = setdiff(1:length(r),rb);
        end
    end
    
end

