classdef SaccadeDetector < handle
    %SACCADEDETECTOR Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Jorge Otero-Millan, jorgeoteromillan@gmail.com 2/17/2014
    %
    
    properties
        MIN_ISI = 20;
        SACCADELIMITTH = 3;
        SACCADELIMITTHACC = 0;
    end
    
    methods(Static=true)
        
        function saccadeDetector = getSaccadeDetector( )
            saccadeDetector = ClusterDetection.SaccadeDetectorCluster( );
            saccadeDetector.Init( );
        end
    end
    
    methods        
        function sac = FindValidSaccades(this, eyeRecording, sac)
            %-- remove sacc in intertrials and in blinks
            valid = eyeRecording.valid;
            
            if ( ~isempty( sac ) )
                
                is_good = zeros(size(sac,1),1);
                for i=1:length(is_good)
                    is_good(i) = sum( ~valid( sac(i,1):sac(i,2)) ) == 0;
                end
                
                sac = sac( is_good==1, :);
                
            end
        end
        
        function sac = FindSaccadeLimits(this, eyeRecording, saccadesIdx)
            % find the begining and the end of the saccade
            
            pi = saccadesIdx;
            starts = zeros(size(pi));
            ends = zeros(size(pi));
            trials = eyeRecording.trials;
            
            v = eyeRecording.getAvgPolarVelocityAcrossEyes();
            acc = eyeRecording.getAvgPolarAccelerationAcrossEyes();
            
            for ipeak=1:length(pi)
                prevtrial = max(trials(trials(:,1)<=pi(ipeak),1));
                nextrial = min(trials(trials(:,2)>=pi(ipeak),2));
                if ( 0 )
                    % this would find the begining only after the
                    % previous peak and the end only until the next
                    % peak
                    if ( ipeak==1)
                        beg = trials(1,1);
                    else
                        beg = max(prevtrial,pi(ipeak-1)+3);
                    end
                    if ( ipeak==length(pi))
                        fin = trials(end,2);
                    else
                        fin = min(nextrial, pi(ipeak+1)-3);
                    end
                else
                    beg = prevtrial;
                    fin = nextrial;
                end
                
                % find the begining
                a = find(v(beg:(pi(ipeak)-1))<this.SACCADELIMITTH | acc(beg:(pi(ipeak)-1))<this.SACCADELIMITTHACC,1,'last');
                if ( ~isempty(a) )
                    a = min(beg + a, pi(ipeak)-1);
                else
                    [mina a] = min(v(beg:(pi(ipeak)-1)));
                    a = beg-1 + a +1;
                end
                starts(ipeak) = a;
                
                
                % find the end
                b = find(v((pi(ipeak)+1):fin)<this.SACCADELIMITTH | acc((pi(ipeak)+1):fin)<this.SACCADELIMITTHACC,1,'first');
                if ( ~isempty(b) )
                    b = max(pi(ipeak)-1 + b, pi(ipeak)+1);
                else
                    [minb b] = min(v((pi(ipeak)+1):fin));
                    b = pi(ipeak) + b -1;
                end

                ends(ipeak) = b;
                
            end
            sac = [starts ends];
            
            %% avoid overlapping peaks, merge the ones that olverlap
            uidx = lohi2idx(sac(:,1), sac(:,2));
            u = zeros(max(sac(:,2))+1,1);
            u(uidx) = 1;
            starts = find(diff([0;u])>0);
            ends = find(diff([u;0])<0);
            
            sac = [starts ends];
        end
        
    end
    
    methods (Static = true)
       
        function sacc_props_enum = GetEnum()
            
            enum.startIndex = 1;
            enum.endIndex = 2;
            enum.duration = 3;
            
            enum.amplitude = 4;
            enum.leftAmplitude = 5;
            enum.rightAmplitude = 6;
            
            enum.displacement = 7;
            enum.leftDisplacement = 8;
            enum.rightDisplacement = 9;
            
            enum.peakVelocity = 10;
            enum.leftPeakVelocity = 11;
            enum.rightPeakVelocity = 12;
            
            enum.peakVelocityIdx = 13;
            enum.leftPeakVelocityIdx = 14;
            enum.rightPeakVelocityIdx = 15;
            
            enum.meanVelocity = 16;
            enum.leftMeanVelocity = 17;
            enum.rightMeanVelocity = 18;
            
            enum.peakAcceleration = 19;
            enum.leftPeakAcceleration = 20;
            enum.rightPeakAcceleration = 21;
            
            enum.peakAccelerationIdx = 22;
            enum.leftPeakAccelerationIdx = 23;
            enum.rightPeakAccelerationIdx = 24;
            
            enum.peakAccelerationStart = 25;
            enum.leftPeakAccelerationStart = 26;
            enum.rightPeakAccelerationStart = 27;
            
            enum.peakAccelerationStartIdx = 28;
            enum.leftPeakAccelerationStartIdx = 29;
            enum.rightPeakAccelerationStartIdx = 30;
            
            enum.peakAccelerationBrake = 31;
            enum.leftPeakAccelerationBrake = 32;
            enum.rightPeakAccelerationBrake = 33;
            
            enum.peakAccelerationBrakeIdx = 34;
            enum.leftPeakAccelerationBrakeIdx = 35;
            enum.rightPeakAccelerationBrakeIdx = 36;

            enum.direction = 37;
            enum.leftDirection = 38;
            enum.rightDirection = 39;
            
            enum.preISI = 40;
            enum.postISI = 41;
         
            sacc_props_enum = enum;
        end
        
        function sacc_props = GetSaccadeProps( eyeRecording, sac  )
            import ClusterDetection.SaccadeDetector;
            
            enum = SaccadeDetector.GetEnum();
            
            sacc_props = zeros( size(sac,1), 40 );
            
            if isempty(sac) 
                return
            end
            
            pos = eyeRecording.position;
            vel = eyeRecording.velocity;
            acc = eyeRecording.acceleration;
            
            posAvg = eyeRecording.getAvgPositionAcrossEyes();
            velAvg = eyeRecording.getAvgPolarVelocityAcrossEyes();
            accAvg = eyeRecording.getAvgPolarAccelerationAcrossEyes();
            
            
            sacc_props(:,enum.startIndex) = sac(:,1);
            sacc_props(:,enum.endIndex) = sac(:,2);
            sacc_props(:,enum.duration) = (sac(:,2) - sac(:,1)) * 1000 / eyeRecording.samplerate;
            
            sacc_props(:,enum.leftAmplitude) = SaccadeDetector.GetAmplitude( pos(:,[eyeRecording.LH eyeRecording.LV]), sac );
            sacc_props(:,enum.rightAmplitude) = SaccadeDetector.GetAmplitude( pos(:,[eyeRecording.RH eyeRecording.RV]), sac );
            sacc_props(:,enum.amplitude) = SaccadeDetector.GetAmplitude( posAvg, sac );
            
            sacc_props(:,enum.leftDisplacement) = SaccadeDetector.GetDisplacement( pos(:,[eyeRecording.LH eyeRecording.LV]), sac );
            sacc_props(:,enum.rightDisplacement) = SaccadeDetector.GetDisplacement( pos(:,[eyeRecording.RH eyeRecording.RV]), sac );
            sacc_props(:,enum.displacement) = SaccadeDetector.GetDisplacement( posAvg, sac );
            
            [sacc_props(:,enum.leftPeakVelocity) sacc_props(:,enum.leftPeakVelocityIdx)] = SaccadeDetector.GetPeakVelocity( vel(:,[eyeRecording.LP]), sac );
            [sacc_props(:,enum.rightPeakVelocity) sacc_props(:,enum.rightPeakVelocityIdx)] = SaccadeDetector.GetPeakVelocity( vel(:,[eyeRecording.RP]), sac );
            [sacc_props(:,enum.peakVelocity) sacc_props(:,enum.peakVelocityIdx)] = SaccadeDetector.GetPeakVelocity(velAvg, sac );
            
            sacc_props(:,enum.leftMeanVelocity) = SaccadeDetector.GetMeanVelocity( vel(:,[eyeRecording.LP]), sac );
            sacc_props(:,enum.rightMeanVelocity) = SaccadeDetector.GetMeanVelocity( vel(:,[eyeRecording.RP]), sac );
            sacc_props(:,enum.meanVelocity) = SaccadeDetector.GetMeanVelocity( velAvg, sac );
            
            [sacc_props(:,enum.leftPeakAcceleration) sacc_props(:,enum.leftPeakAccelerationIdx)] = SaccadeDetector.GetPeakAcceleration( acc(:,[eyeRecording.LP]), sac);
            [sacc_props(:,enum.rightPeakAcceleration) sacc_props(:,enum.rightPeakAccelerationIdx)] = SaccadeDetector.GetPeakAcceleration( acc(:,[eyeRecording.RP]), sac);
            [sacc_props(:,enum.peakAcceleration) sacc_props(:,enum.peakAccelerationIdx)] = SaccadeDetector.GetPeakAcceleration( accAvg, sac );
            
            [sacc_props(:,enum.leftPeakAccelerationStart) sacc_props(:,enum.leftPeakAccelerationStartIdx)] = SaccadeDetector.GetPeakAccelerationStart( acc(:,[eyeRecording.LP]), sac, sacc_props(:,enum.leftPeakVelocityIdx));
            [sacc_props(:,enum.rightPeakAccelerationStart) sacc_props(:,enum.rightPeakAccelerationStartIdx)] = SaccadeDetector.GetPeakAccelerationStart( acc(:,[eyeRecording.RP]), sac, sacc_props(:,enum.rightPeakVelocityIdx));
            [sacc_props(:,enum.peakAccelerationStart) sacc_props(:,enum.peakAccelerationStartIdx)] = SaccadeDetector.GetPeakAccelerationStart( accAvg, sac, sacc_props(:,enum.peakVelocityIdx));
            
            [sacc_props(:,enum.leftPeakAccelerationBrake) sacc_props(:,enum.leftPeakAccelerationBrakeIdx)] = SaccadeDetector.GetPeakAccelerationBrake( acc(:,[eyeRecording.LP]), sac, sacc_props(:,enum.leftPeakVelocityIdx));
            [sacc_props(:,enum.rightPeakAccelerationBrake) sacc_props(:,enum.rightPeakAccelerationBrakeIdx)] = SaccadeDetector.GetPeakAccelerationBrake( acc(:,[eyeRecording.RP]), sac, sacc_props(:,enum.rightPeakVelocityIdx));
            [sacc_props(:,enum.peakAccelerationBrake) sacc_props(:,enum.peakAccelerationBrakeIdx)] = SaccadeDetector.GetPeakAccelerationBrake( accAvg, sac, sacc_props(:,enum.peakVelocityIdx));

            sacc_props(:,enum.leftDirection) = SaccadeDetector.GetDirection( pos(:,[eyeRecording.LH eyeRecording.LV]), sac );
            sacc_props(:,enum.rightDirection) = SaccadeDetector.GetDirection( pos(:,[eyeRecording.LH eyeRecording.RV]), sac );
            sacc_props(:,enum.direction) = SaccadeDetector.GetDirection( posAvg, sac );
            
            
            [sacc_props(:,enum.preISI) sacc_props(:,enum.postISI)] = SaccadeDetector.GetISI( sac, eyeRecording.valid );
            
        end
        
        function amplitude = GetAmplitude( pos, sac  )
            amplitude = zeros(size(sac,1),1);
            for i = 1:size(sac,1)
                idx = sac(i,1):sac(i,2);
                dX = max(pos(idx,1)) - min(pos(idx,1));
                dY = max(pos(idx,2)) - min(pos(idx,2));
                amplitude(i)  = sqrt( dX.^2 + dY.^2);
            end
        end
        
        function displacement = GetDisplacement( pos, sac  )
            displacement = sqrt( (pos(sac(:,2),1) - pos(sac(:,1),1) ).^2 + ( pos(sac(:,2),2) - pos(sac(:,1),2) ).^2 );
        end
        
        function [pkvel idxmax] = GetPeakVelocity( vel, sac  )
            pkvel = zeros(size(sac,1),1);
            idxmax = zeros(size(sac,1),1);
            for i = 1:size(sac,1)
                idx = sac(i,1):sac(i,2);
                [pkvel(i) idxmax(i)] = max( vel(idx) );
            end
        end
        
        function mnvel = GetMeanVelocity( vel, sac  )
            mnvel = zeros(size(sac,1),1);
            for i = 1:size(sac,1)
                idx = sac(i,1):sac(i,2);
                mnvel(i) = mean( vel(idx) );
            end
        end
        
        function [pkacc idxmax] = GetPeakAcceleration( acc, sac  )
            pkacc = zeros(size(sac,1),1);
            idxmax = zeros(size(sac,1),1);
            for i = 1:size(sac,1)
                idx = sac(i,1):sac(i,2);
                [pkacc(i) idxmax(i)] = max( acc(idx) );
            end
        end
        
        function  [pkacc idxmax] = GetPeakAccelerationStart( acc, sac, sacPeakVelocityIdx  )
            pkacc = zeros(size(sac,1),1);
            idxmax = zeros(size(sac,1),1);
            for i = 1:size(sac,1)
                idx = sac(i,1):(sac(i,1)+sacPeakVelocityIdx(i));
                [pkacc(i) idxmax(i)] = max( acc(idx) );
            end
        end
        
        function [pkacc idxmax]  = GetPeakAccelerationBrake( acc, sac, sacPeakVelocityIdx  )
            pkacc = zeros(size(sac,1),1);
            idxmax = zeros(size(sac,1),1);
            for i = 1:size(sac,1)
                idx = min((sac(i,1)+sacPeakVelocityIdx(i)),sac(i,2)):sac(i,2);
                [pkacc(i) idxmax(i)] = max( acc(idx) );
                idxmax(i) = idxmax(i) + idx(1)-1;
            end
        end
        
        function direction = GetDirection( pos, sac )
            dx = pos(sac(:,2),1) - pos(sac(:,1),1);
            dy = pos(sac(:,2),2) - pos(sac(:,1),2);
            directions	= atan2( dx, dy );
            direction = mod( (directions * (180 / pi)) + 360, 360);
        end
        
        function [preisi postisi] = GetISI( sac, valid )
            if ( isempty( sac) )
                preisi = zeros(0,1);
                postisi = zeros(0,1);
                return;
            end
            preisi = nan(size(sac,1),1);
            postisi = nan(size(sac,1),1);
            for i = 2:size(sac,1)
                currbegin = sac(i,1);
                lastend = max(sac(sac(:,2)<=currbegin,2));
                if ( sum(~valid(lastend:currbegin)) > 0 )
                    preisi(i) = nan;
                    postisi(i-1) = nan;
                else
                    preisi(i) = currbegin-lastend+1;
                    postisi(i-1) = currbegin-lastend+1;
                end
            end
            
        end
    end
    
end

