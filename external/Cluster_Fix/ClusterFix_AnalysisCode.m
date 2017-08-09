%% Code runs and plots several key components of Cluster Fix Algorithm
% Written by Seth Konig 2013. Buffalo Lab.
%
% Searching [#] will bring you to that # section
%
% 1. Run Cluster Fix
% 2. Plot Scan Paths with fixations in red and saccades in green
% 3. Plot Velocity Traces with fixations in red and saccades in green
% 4. Plot Acceleration Traces with fixations in red and saccades in green
% 5. Plot Global Clusters in 3D state space by global cluster
% 6. Plot Scan Path by global cluster
% 7. Plot Velocity and Acceleration by global cluster
% 8. Plot Local Clusters in 3D state space by global cluster
% 9. Plot Scan Path by local cluster
% 10. Plot Velocity and Acceleration by local cluster
% 11. Plot Scan Path after combing local re-clustering with global results
% 12. Plot Velocity and Acceleration after combing local re-clustering with global results
% 13. Plot Global State Space after combing local re-clustering with global results
% 14. Support Vector Machine: analyze seperatibility of fixations from saccades
%
% Files you need:
% 1. MP100712_1.mat: raw eye traces that occured within image boundary
%    Contains 36 eye traces from the viewing of 36 novel image.
% 2. ClusterFix.m: Cluster Fix algorithm
%
% All units are in dva (degrees of visual angle)
%% [1] Run Cluster Fix

load('MP100712_1','eyedat') %load eye data cell array
fixationstats = ClusterFix(eyedat);

%% [2] Plot Scan Paths with fixations in red and saccades in green
% scale in degrees of visual angle (dva)

for i = 1:length(fixationstats);
    xy = fixationstats{i}.XY;
    fixations = fixationstats{i}.fixations;
    fixationtimes = fixationstats{i}.fixationtimes;
    
    figure
    hold on
    plot(xy(1,:),xy(2,:),'g'); %green for saccades
    for ii = 1:length(fixationtimes);
        plot(xy(1,fixationtimes(1,ii):fixationtimes(2,ii)),...
            xy(2,fixationtimes(1,ii):fixationtimes(2,ii)),'r'); %red for fixations
        plot(fixations(1,ii),fixations(2,ii),'k*'); %plot mean fixation location
    end
    legend('Saccades','Fixations')
    
end
%% [3] Plot Velocity Traces with fixations in red and saccades in green

fltord = 60; %filter order
lowpasfrq = 30; %Low pass frequency cutoff
nyqfrq = 200 ./ 2; %nyquist frequency
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
%note when filtering scan paths for fixation detection we also use an
%upsampling process that helps increase the signal-to-noise ratio

for i = 1:length(fixationstats);
    fixationtimes = fixationstats{i}.fixationtimes;
    xy = fixationstats{i}.XY;
    
    x = xy(1,:);
    y = xy(2,:);
    x = [x(20:-1:1) x x(end:-1:end-20)]; %add 20 ms buffer for filtering
    y = [y(20:-1:1) y y(end:-1:end-20)]; %add 20 ms buffer for filtering
    x = filtfilt(flt,1,x); %filter
    y = filtfilt(flt,1,y); %filter
    x = x(21:end-20); %remove buffer after filtering
    y = y(21:end-20); %remove buffer after filtering
    
    velx = diff(xy(1,:)); %x-component of velocity
    vely = diff(xy(2,:)); %y-component of velocity
    vel = sqrt(velx.^2 + vely.^2);
    vel = vel*200;%convert to dva  per sec since sampled at 200 Hz
    
    figure
    plot(vel,'g'); %green for saccades
    hold on
    for ii = 1:size(fixationtimes,2);
        fixtimes = fixationtimes(1,ii):fixationtimes(2,ii);
        fixtimes(fixtimes > length(vel)) = length(vel);
        plot(fixtimes,vel(fixtimes),'r'); %red for fixations
    end
    legend('Saccades','Fixations','location','NorthEastOutside')
    xlabel('Sample (5 ms/sample)')
    ylabel('Velocity (dva/sec)')
    box off
    xlim([1 500])
end
%% [4] Plot Acceleration Traces with fixations in red and saccades in green

fltord = 60; %filter order
lowpasfrq = 30; %Low pass frequency cutoff
nyqfrq = 200 ./ 2; %nyquist frequency
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
%note when filtering scan paths for fixation detection we also use an
%upsampling process that helps increase the signal-to-noise ratio

for i = 1:length(fixationstats);
    fixationtimes = fixationstats{i}.fixationtimes;
    xy = fixationstats{i}.XY;
    
    x = xy(1,:);
    y = xy(2,:);
    x = [x(20:-1:1) x x(end:-1:end-20)]; %add 20 ms buffer for filtering
    y = [y(20:-1:1) y y(end:-1:end-20)]; %add 20 ms buffer for filtering
    x = filtfilt(flt,1,x); %filter
    y = filtfilt(flt,1,y); %filter
    x = x(21:end-20); %remove buffer after filtering
    y = y(21:end-20); %remove buffer after filtering
    
    velx = diff(xy(1,:));
    vely = diff(xy(2,:));
    vel = sqrt(velx.^2 + vely.^2);
    vel = vel;%convert to dva  per sec since sampled at 200 Hz
    accel = abs(diff(vel))*200*200; %in dva/sec/sec
    
    figure
    plot(accel,'g');
    hold on
    for ii = 1:size(fixationtimes,2);
        fixtimes = fixationtimes(1,ii):fixationtimes(2,ii);
        fixtimes(fixtimes > length(accel)) = length(accel);
        plot(fixtimes,accel(fixtimes),'r'); %red for fixations
    end
    legend('Saccades','Fixations','location','NorthEastOutside')
    xlabel('Sample (5 ms/sample)')
    ylabel('Acceleration (dva/sec^2)')
    box off
    xlim([1 500])
end
%% [5] Plot Global Clusters in 3D state space by global cluster
% Plots global clustering results prior to local re-clustering
% Run ClusterFix in debug mode with a break point at line 110 or line 131
% Run to line 110: Clusters determined by k-means clustering
% Run to line 131: Clusters consolidated into fixations (red) vs saccades (green)

%uses weird variable names because in debug mode and to not overwrite important variables
ia = 'rgbmk';
figure
hold on
for TIND = 1:max(T);
    plot3(points((T == TIND),2),points((T == TIND),3),points((T == TIND),4),...
        [ia(TIND)  '.'],'markersize',6); %points((T == TIND),1) will plot distance
end
xlabel('velocity')
ylabel('acceleration')
zlabel('Angular Velocity')
view(-7,32)
%% [6] Plot Scan Path by global cluster
% Plots scan path (x's and y's) for global clustering results prior to local re-clustering
% Run ClusterFix in debug mode with a break point at line 110 or line 131
% Run to line 110: Clusters determined by k-means clustering
% Run to line 131: Clusters consolidated into fixations (red) vs saccades (green)

%uses weird variable names because in debug mode and to not overwrite important variables
x = eyedat{cndlop}(1,:);
y = eyedat{cndlop}(2,:);
x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
y = [y(buffer:-1:1) y y(end:-1:end-buffer)];  %add buffer for filtering
xss = resample(x,samprate*1000,1);%up sample to 1000 Hz
yss = resample(y,samprate*1000,1);%up sample to 1000 Hz
xss = xss(101:end-100); %remove buffer
yss = yss(101:end-100); %remove buffer

a = 1:length(xss);
ia = 'rgbmk';
figure
hold on
plot(xss,yss,'k')
for TIND = 1:max(T);
    x = a(T == TIND);
    y = find(diff(x) > 1);
    y =  [1 y length(x)];
    for ii = 2:length(y);
        rng = x(y(ii-1)+1:y(ii));
        plot(xss(rng),yss(rng),ia(TIND))
    end
end

% Should re-run the following afterwards so you don't mess up these variables
x = eyedat{cndlop}(1,:);
y = eyedat{cndlop}(2,:);
x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
x = resample(x,samprate*1000,1);%up sample to 1000 Hz
y = resample(y,samprate*1000,1);%up sample to 1000 Hz
xss = filtfilt(flt,1,x);
yss = filtfilt(flt,1,y);
xss = xss(101:end-100); %remove buffer after filtering
yss = yss(101:end-100); %remove buffer after filtering
x = x(101:end-100); %remove buffer after filtering
y = y(101:end-100); %remove buffer after filtering
%% [7] Plot Velocity and Acceleration by global cluster
% Plots velocity and acceleration profiles for global clustering results prior to local re-clustering
% Run ClusterFix in debug mode with a break point at line 110 or line 131
% Run to line 110: Clusters determined by k-means clustering
% Run to line 131: Clusters consolidated into fixations (red) vs saccades (green)

%uses weird variable names because in debug mode and to not overwrite important variables
a = 1:length(xss);
ia = 'rgbmk';
figure
hold on
plot(vel,'k')
for TIND = 1:max(T);
    x = a(T == TIND);
    y = find(diff(x) > 1);
    y =  [1 y length(x)];
    for ii = 2:length(y);
        rng = x(y(ii-1)+1:y(ii));
        plot(rng,vel(rng),ia(TIND))
    end
end
xlabel('Time (ms)')
ylabel('Normalized Velocity')

a = 1:length(xss);
ia = 'rgbmk';
figure
hold on
plot(accel,'k') %background because points are disjointed
for TIND = 1:max(T);
    x = a(T == TIND);
    y = find(diff(x) > 1);
    y =  [1 y length(x)];
    for ii = 2:length(y);
        rng = x(y(ii-1)+1:y(ii));
        plot(rng,accel(rng),ia(TIND))
    end
end
xlabel('Time (ms)')
ylabel('Normalized Acceleration')

% Should re-run the following afterwards so you don't mess up these variables
x = eyedat{cndlop}(1,:);
y = eyedat{cndlop}(2,:);
x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
x = resample(x,samprate*1000,1);%up sample to 1000 Hz
y = resample(y,samprate*1000,1);%up sample to 1000 Hz
xss = filtfilt(flt,1,x);
yss = filtfilt(flt,1,y);
xss = xss(101:end-100); %remove buffer after filtering
yss = yss(101:end-100); %remove buffer after filtering
x = x(101:end-100); %remove buffer after filtering
y = y(101:end-100); %remove buffer after filtering
%%  [8] Plot Local Clusters in 3D state space by local cluster
% Plots local clustering results after initial global clustering
% Run ClusterFix in debug mode with a break point at line 152 or line 179
% Run to line 152: Clusters determined by k-means clustering
% Run to line 179: Clusters consolidated into fixations (red) vs saccades (green)
%
% Since this is local re-evalution step can run this FOR loop (starting at
% line 137 and ending at 180) for each locally re-evaluated fixation

%uses weird variable names because in debug mode and to not overwrite important variables

ia = 'rgbmk';
figure
hold on
for iii = 1:max(T);
    pp = find(T == iii);
    plot3(POINTS(pp,2),POINTS(pp,3),POINTS(pp,4),[ia(iii)  '.']); %POINTS(pp,1) will plot distance
end
view(3)
xlabel('velocity')
ylabel('acceleration')
zlabel('Angular Velocity')
%% [9] Plot Scan Path by local cluster
% Plots local clustering results after initial global clustering
% Run ClusterFix in debug mode with a break point at line 152 or line 179
% Run to line 152: Clusters determined by k-means clustering
% Run to line 179: Clusters consolidated into fixations (red) vs saccades (green)
%
% Since this is local re-evalution step can run this FOR loop (starting at
% line 137 and ending at 180) for each locally re-evaluated fixation

% uses weird variable names because in debug mode and to not overwrite important variables

ia = 'rgbmk';
figure
hold on
plot(x(altind),y(altind),'k') %background because points are disjointed
for ii = 1:max(T);
    pp = altind(find(T == ii));
    ib = find(diff(pp) > 1);
    if isempty(ib);
        plot(x(pp),y(pp),ia(ii));
    else
        for iii = 1:length(ib)+1;
            if iii == 1;
                plot(x(pp(1:ib(iii))),y(pp(1:ib(iii))),ia(ii));
            elseif iii == length(ib)+1;
                plot(x(pp(ib(end)+1:end)),y(pp(ib(end)+1:end)),ia(ii));
            else
                plot(x(pp(ib(iii-1)+1:ib(iii))),y(pp(ib(iii-1)+1:ib(iii))),ia(ii));
            end
        end
    end
end
%% [10] Plot Velocity and Acceleration by local cluster
% Plots local clustering results after initial global clustering
% Run ClusterFix in debug mode with a break point at line 152 or line 179
% Run to line 152: Clusters determined by k-means clustering
% Run to line 179: Clusters consolidated into fixations (red) vs saccades (green)
%
% Since this is local re-evalution step can run this FOR loop (starting at
% line 137 and ending at 180) for each locally re-evaluated fixation

% uses weird variable names because in debug mode and to not overwrite important variables

% points velocity is in column 2 and acceleration in column 3 of variable points
vel = vel/(mean(vel)+3*std(vel));
ia = 'rgbmk';
figure
hold on
plot(altind,vel(altind),'k') %background because points are disjointed
for ii = 1:max(T);
    pp = altind(find(T == ii));
    ib = find(diff(pp) > 1);
    if isempty(ib);
        plot(pp,vel(pp),ia(ii));
    else
        for iii = 1:length(ib)+1;
            if iii == 1;
                plot(pp(1:ib(iii)),vel(pp(1:ib(iii))),ia(ii));
            elseif iii == length(ib)+1;
                plot(pp(ib(end)+1:end),vel(pp(ib(end)+1:end)),ia(ii));
            else
                plot(pp(ib(iii-1)+1:ib(iii)),vel(pp(ib(iii-1)+1:ib(iii))),ia(ii));
            end
        end
    end
end
xlabel('Time (ms)')
ylabel('Normalized Velocity')

accel = accel/(mean(accel)+3*std(accel));
ia = 'rgbmk';
figure
hold on
plot(altind,accel(altind),'k') %background because points are disjointed
for ii = 1:max(T);
    pp = altind(find(T == ii));
    ib = find(diff(pp) > 1);
    if isempty(ib);
        plot(pp,accel(pp),ia(ii));
    else
        for iii = 1:length(ib)+1;
            if iii == 1;
                plot(pp(1:ib(iii)),accel(pp(1:ib(iii))),ia(ii));
            elseif iii == length(ib)+1;
                plot(pp(ib(end)+1:end),accel(pp(ib(end)+1:end)),ia(ii));
            else
                plot(pp(ib(iii-1)+1:ib(iii)),accel(pp(ib(iii-1)+1:ib(iii))),ia(ii));
            end
        end
    end
end
xlabel('Time (ms)')
ylabel('Normalized Acceleration')
%% [11] Plot Scan Path after combing local re-clustering with global results
% Plots global clusters after local-reclustering with fixations (red) vs saccades (green)
% Run ClusterFix in debug mode with a break point at line 219

figure
hold on
for ii = 1:size(saccadetimes,2);
    plot(x(saccadetimes(1,ii):saccadetimes(2,ii)),...
        y(saccadetimes(1,ii):saccadetimes(2,ii)),'g')
end
for ii = 1:size(fixationtimes,2);
    plot(x(fixationtimes(1,ii):fixationtimes(2,ii)),...
        y(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
end
hold off

%% [12] Plot Velocity and Acceleration after combing local re-clustering with global results
% Plots global clusters after local-reclustering with fixations (red) vs saccades (green)
% Run ClusterFix in debug mode with a break point at line 219

% uses weird variable names because in debug mode and to not overwrite important variables

figure
hold on
for ii = 1:size(saccadetimes,2);
    plot(saccadetimes(1,ii):saccadetimes(2,ii),...
        1000*vel(saccadetimes(1,ii):saccadetimes(2,ii)),'g') %*1000 to conver to dva/sec
end
for ii = 1:size(fixationtimes,2);
    plot(fixationtimes(1,ii):fixationtimes(2,ii),...
        1000*vel(fixationtimes(1,ii):fixationtimes(2,ii)),'r') %*1000 to conver to dva/sec
end
hold off
xlabel('Time (ms)')
ylabel('Velocity (dva/sec)')


figure
hold on
for ii = 1:size(saccadetimes,2);
    plot(saccadetimes(1,ii):saccadetimes(2,ii),...
        1000^2*accel(saccadetimes(1,ii):saccadetimes(2,ii)),'g') %*1000^2 to conver to dva/sec/sec
end
for ii = 1:size(fixationtimes,2);
    plot(fixationtimes(1,ii):fixationtimes(2,ii),...
        1000^2*accel(fixationtimes(1,ii):fixationtimes(2,ii)),'r') %*1000^2 to conver to dva/sec/sec
end
hold off
xlabel('Time (ms)')
ylabel('Acceleration (dva/sec^2)')
%% [13] Plot Global State Space after combing local re-clustering with global results
% Plots global clusters after local-reclustering with fixations (red) vs saccades (green)
% Run ClusterFix in debug mode with a break point at line 219
% Also plots velocity and acceleration thresholds (mean + std)

% uses weird variable names because in debug mode and to not overwrite important variables

T = zeros(1,length(points));
T(fixationindexes)=1;
T(T==0) = 2;

ia = 'rg';
figure
hold on
for TIND = 1:max(T);
    plot3(points((T == TIND),2),points((T == TIND),3),points((T == TIND),4),...
        [ia(TIND)  '.'],'markersize',6);
end
xlabel('velocity')
ylabel('acceleration')
zlabel('Angular Velocity')
view(2) % alternatively can use view(-7,32)

plot([0 1],[(mean(accel)+std(accel))/(mean(accel)+3*std(accel))...
    (mean(accel)+std(accel))/(mean(accel)+3*std(accel))],'k'); %acceleration threshold
plot([(mean(vel)+std(vel))/(mean(vel)+3*std(vel))...
    (mean(vel)+std(vel))/(mean(vel)+3*std(vel))],[0 1],'k'); %velocity threshold
%% [14] Support Vector Machine: analyze seperatibility of fixations from saccades
% calculate parameter (function on line 343) values for fixations vs saccades by image
% Must run [1] to get fixationstats.
% MeanClusterValues: row 1 is fixations and row 2 is saccades ;col 1: max vel,
% col 2: max accel, col 3: mean distance, col 4: mean vel, col 5: |mean anlge|, col 6: mean Angular Velocity

fixclusters = [];
sacclusters = [];
for cndlop = 1:length(fixationstats)
    fixclusters = [fixclusters ...
        fixationstats{cndlop}.MeanClusterValues(1,2:2:end)'];
    sacclusters = [sacclusters ...
        fixationstats{cndlop}.MeanClusterValues(2,2:2:end)'];
end

figure
hold on
for ii = 1:size(fixclusters,2)
    plot3(fixclusters(1,ii),fixclusters(2,ii),200*fixclusters(3,ii),'r.')
end
for ii = 1:size(sacclusters,2)
    plot3(sacclusters(1,ii),sacclusters(2,ii),200*sacclusters(3,ii),'g.')
    %*200 to convert from degrees/sample to degrees/sec
end
view(3)
xlabel('Max Accleration (dva/sec^{2})')
ylabel('Velocity (dva/sec)')
zlabel('Angular Velocity (degrees/sec)')
view(65,40)

% measure seperability with a support vector machine (SVM)
correct = NaN(1,100);
totalpoints = NaN(1,100);
for iter = 1:100;
    data = [fixclusters'; sacclusters'];
    groups = [true(1,size(fixclusters,2)) false(1,size(sacclusters,2))]';
    P = cvpartition(groups,'Holdout',0.75); % 0.75 uses 25% of points to train and 75% to test
    SVMstruct = svmtrain(data(P.training,:),groups(P.training));%,'showplot','true');
    testresults = svmclassify(SVMstruct,data(P.test,:));%,'showplot','true');
    correct(iter) = sum(groups(P.test) == testresults);
    totalpoints(iter) = sum(P.test);
end
lowestaccuracy = 100*min(correct./totalpoints)
averageaccuracy = 100*mean(correct./totalpoints)