classdef SaccadeDetectorCluster < ClusterDetection.SaccadeDetector
    %SACCADEDETECTORCLUSTER Implementation of the saccade detection method
    %using kmeans
    %
    % Jorge Otero-Millan, jorgeoteromillan@gmail.com 2/17/2014
    %
    
    properties
        MINIPI = 0; % minimum inter peak interval
        RATEOFPEAKS = 5; % number of peaks per second
        MINPEAKVEL = 1;
        
        featureSelection
        
        NumMaxClusters
    end
    
    methods
        function Init( this )
            this.NumMaxClusters = 4;
            this.featureSelection = {'logmeanvel', 'logmeanaccelerationStart', 'logmeanaccelerationBreak'};
            this.RATEOFPEAKS = 5;
        end
        
        function [sac stats] = FindSaccades( this, eyeRecording ) 
           
            enum = ClusterDetection.SaccadeDetector.GetEnum();
            
            peaks = this.FindPeaks( eyeRecording );
            peaks = this.GetSaccadeProps( eyeRecording, peaks  );

            % get the amplitude of each peak (it will be useful to sort the
            % clusters by average amplitude)
            if ( eyeRecording.hasLeftEye && eyeRecording.hasRightEye )
                vels = peaks(:,enum.peakVelocity);
            elseif( eyeRecording.hasLeftEye  )
                vels = peaks(:,enum.leftPeakVelocity);
            elseif ( eyeRecording.hasRightEye )
                vels = peaks(:,enum.rightPeakVelocity);
            end
            
            [features featuresBeforeWhiten featuresRotated allfeatures] = this.GetFeatures( eyeRecording, peaks  );
            stats.features = features;
            stats.featuresBeforeWhiten = featuresBeforeWhiten;
            stats.featuresRotated = featuresRotated;
            stats.allfeatures = allfeatures;
            stats.featureselection = this.featureSelection;
            stats.peaks = peaks;
           
            trials= eyeRecording.trials;
            
            % we will cluster trial by trial unless the trial has less than
            % 100 events, in that case we will group trials until we have
            % more than 100. 
            % I think this is specially good just in case the properties of
            % the recording change over time
            clusteridx = [];
            chunkNumber = [];
            currentChunk = 1;
            ctrs = {};
            iTrial = 0;
            while iTrial < length(trials)
                idx = [];
                while(length(idx) < 500 )
                    
                    iTrial = iTrial + 1;
                    idxInTrial = find ( trials(iTrial,1) <= peaks(:,enum.startIndex) & peaks(:,enum.startIndex) < trials(iTrial,2));
                    idx = cat(1,idx, idxInTrial);
                    
                    % this is to take care of the last chunk of data having
                    % less than 100 events. If that's the case we take them
                    % now
                    idxLeft = find( peaks(:,enum.startIndex) >= trials(iTrial,2) );
                    if ( length(idxLeft) < 100 )
                        idx = cat(1,idx, idxLeft);
                        iTrial = length(trials) + 1;
                        break
                    end
                end
                
                
                [clusteridxTrial ctrsTrial] = this.ClusterPeaks( features(idx,:), vels(idx) );
                
                clusteridx = cat(1, clusteridx, clusteridxTrial);
                
                chunkNumber = cat(1, chunkNumber, currentChunk*ones(size(clusteridxTrial)));
                currentChunk = currentChunk+1;
                
                ctrs{end+1} = ctrsTrial;

            end
            
            stats.silhouette = mean(silhouette(features, min(clusteridx,2)));
            stats.clusteridx = clusteridx;
            stats.ctrs = ctrs;
            stats.chunkNumber = chunkNumber;
            
            sac = peaks(clusteridx==1, :);
            sac = this.GetSaccadeProps( eyeRecording, sac  );
        end
        
        function peaks = FindPeaks( this, eyeRecording )
            
            this.MINIPI = 20*eyeRecording.samplerate/1000; % minimum inter peak interval
                        
            vel = eyeRecording.getAvgPolarVelocityAcrossEyes();
            valid = eyeRecording.valid;

            peaksIdx = [];
            % find the peaks in each trial
            
            trials= eyeRecording.trials;
            
            for iTrial = 1:length(trials)

                idxtrial = trials(iTrial,1):trials(iTrial,2);
                nsamples = length(idxtrial);
                
                % determine the number of peaks for this chunk
                npeaks = ceil((nsamples-sum(~valid(idxtrial)))/eyeRecording.samplerate*this.RATEOFPEAKS);
                if ( npeaks ==0)
                    continue
                end
                
                pidx = this.myfindpeaks(vel(idxtrial)', this.MINIPI,npeaks);

                peaksIdx = cat(1, peaksIdx, pidx' + trials(iTrial,1)-1);
            end
            
            peaks = this.FindSaccadeLimits( eyeRecording, peaksIdx);
            peaks = this.FindValidSaccades( eyeRecording, peaks );
        end
        
        function  peakidx = myfindpeaks( this, data, mindistance, numpeaks )
            peakvalues = zeros(1,numpeaks);
            peakidx = zeros(1,numpeaks);
            
            % find all the local maxima
            allpeaks = find( diff(data(1:end-1))>=0 & diff(data(2:end))<=0 )+1;
            allpeaksvalues = data(allpeaks);
            i=1;
            while(i<=numpeaks)
                % get the highest peak left
                [p pi] = max(allpeaksvalues);
                
                % get the neighborhood data
                endL = max(1,allpeaks(pi) - mindistance);
                endR = min(allpeaks(pi) + mindistance,length(data));
                
                % if the peak is not the highest point within the
                % neighborhood discard it
                if ( ~(any(p < data(endL:endR)) ))
                    peakidx(i) = allpeaks(pi);
                    peakvalues(i) = allpeaksvalues(pi);
                    % if the peak is good discard other peaks within the
                    % neigborhood
                    remidx = find(abs(allpeaks-allpeaks(pi)) <= mindistance);
                    i=i+1;
                else
                    remidx = pi;
                end
                
                allpeaks(remidx) = [];
                allpeaksvalues(remidx) = [];
            end
            % return peaks and peak values sorted
            peakidx = sort(peakidx,2,'ascend');
        end
        
        function [features featuresBeforeWhiten featuresRotated allfeatures ] = GetFeatures( this, eyeRecording, sacc_props )
            
            enum = ClusterDetection.SaccadeDetector.GetEnum();
            
            if ( eyeRecording.hasLeftEye && eyeRecording.hasRightEye )
                features.logmeanvel = log(sacc_props(:,enum.peakVelocity));
                features.logmeanmag = log(sacc_props(:,enum.amplitude));
                features.meanacceleration = log(sacc_props(:,enum.peakAcceleration));
                features.logmeanaccelerationStart = log(sacc_props(:,enum.peakAccelerationStart));
                features.logmeanaccelerationBreak = log(sacc_props(:,enum.peakAccelerationBrake));
            elseif( eyeRecording.hasLeftEye  )
                features.logmeanvel = log(sacc_props(:,enum.leftPeakVelocity));
                features.logmeanmag = log(sacc_props(:,enum.leftAmplitude));
                features.meanacceleration = log(sacc_props(:,enum.leftPeakAcceleration));
                features.logmeanaccelerationStart = log(sacc_props(:,enum.peakAccelerationStart));
                features.logmeanaccelerationBreak = log(sacc_props(:,enum.peakAccelerationBrake));         
            elseif ( eyeRecording.hasRightEye )
                features.logmeanvel = log(sacc_props(:,enum.rightPeakVelocity));
                features.logmeanmag = log(sacc_props(:,enum.rightAmplitude));
                features.meanacceleration = log(sacc_props(:,enum.rightPeakAcceleration));
                features.logmeanaccelerationStart = log(sacc_props(:,enum.peakAccelerationStart));
                features.logmeanaccelerationBreak = log(sacc_props(:,enum.peakAccelerationBrake));
            end
            
            allfeatures = features;
            
            nfeatures = length(this.featureSelection);
            
            % build feature matrix
            X = [];
            for i=1:nfeatures
                f = features.(this.featureSelection{i});
                f(f==-Inf) = min(f(f~=-Inf));
            
                X = cat(2,X, zscore(f));
            end
            
            featuresBeforeWhiten = X;
            
            C= cov(X);
            M= mean(X);
            [V,D]= eig(C);
            P = V * diag(sqrt(1./(diag(D) + 0.1))) ;
            W = bsxfun(@minus, X, M) * P;
            featuresRotated = W;
            
            d = diag(D)/D(end);
            X = W(:,d>0.05);
            features = X;
        end
        
        function [clusteridx ctrs] = ClusterPeaks( this, features, vels )
            
            
            X  = features;
            
            smax = 0;
            % detect 2 clusters, then 3, then 4 ... up to Nclusters
            for Nc=min(2,this.NumMaxClusters):this.NumMaxClusters
                
                opts = statset('Display','off');
                
                starts = zeros(Nc,size(X,2));
                [velssorted velsrank] = sort(vels);
                for i=1:Nc
                    idxstart = (1:floor(length(vels)/Nc)) + floor(length(vels)/Nc)*(i-1);
                    starts(i,:) = mean(X(velsrank(idxstart),:));
                end
                
                [idxtemp,ctrstemp] = kmeans(X(1:end,:),Nc, 'Distance','sqeuclid', 'Options',opts, 'Start',starts);

                Nclusters = max(idxtemp);
                
                % calculate the mean magnitude of each cluster
                Mcluster = zeros(1,Nclusters);
                for i=1:Nclusters
                    Mcluster(i) =  mean(vels((idxtemp==i)));
                end
                
                % sort clusters by mean magnitude
                [MclusterSorted MclusterRank] = sort(Mcluster,2, 'descend');
                idx_temp2 = idxtemp;
                
                % renumber the clusters ordered bhy magnitude,  1 is the
                % one with largest magnitude (the one with the saccades)
                for i=1:Nclusters
                    idxtemp(idx_temp2==MclusterRank(i)) = i;
                end
                
                
                % calculate the mean silhouette but just with 2 clusters,
                % saccades and no-saccades. Grouping all the no-saccades
                % into one cluster with index 2
                s = mean(silhouette(X, min(idxtemp,2)));
                if ( smax*1.01 > s )
                    % if for the current Nc the mean silhouette is smaller
                    % than for the previous Nc break and keep the previous
                    % Nc (number of clusters)
                    break;
                end
                smax = s;
                idx = idxtemp;
                ctrs = ctrstemp;
            end
            
            clusteridx = idx;
            
        end
        
        
        
        
    end
    
end