%STANDARD PROGRAM FOR DATA ANALYSIS: analyze the file "spikes.txt"

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

i2mPath;
BinTab=[10 20];%Put here the different bin sizes you want to test
SizeMax=5;
nbpatterns=1000;
n_sample=1000;
N=10;%Enter here the number of neurons to be included in the analysis
PrintFig=1;%Put to 1 to store some figures in the Figures directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initializing storage matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

divjs=zeros(2,length(BinTab),SizeMax,3);
divjsNoJ1=zeros(2,length(BinTab),SizeMax,3);
divjsNoJ=zeros(2,length(BinTab),SizeMax,3);

coefsh= zeros(N,3,length(BinTab));
coefsJ = zeros(N,N,3,length(BinTab));
coefsJ1 = zeros(N,N,3,length(BinTab));
coefshr = zeros(N,3,length(BinTab));
coefsJr = zeros(N,N,3,length(BinTab));
coefshstat = zeros(N,3,length(BinTab));
coefsJstat = zeros(N,N,3,length(BinTab));

PexpTab = zeros(nbpatterns,length(BinTab),SizeMax,3);
PthTab = zeros(nbpatterns,length(BinTab),SizeMax,3);
PthNoJ1Tab = zeros(nbpatterns,length(BinTab),SizeMax,3);
PthNoJTab = zeros(nbpatterns,length(BinTab),SizeMax,3);

r_raw = LoadRaster('spikes.txt');%Enter here the file where the spike times are stored
f=find(r_raw(2,:)<=N);%Reduce the number of neurons to be included in the analysis to N
r = r_raw(:,f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop on the time bin size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for BinIndex=1:length(BinTab)
    BinSize=BinTab(BinIndex);
    
    % Load the raster data: it is also possible to concatenate the
    % information from several rasters, and to have different types of
    % data (indexed by the variable "state"). 
    raster{1,1} = BinRaster(r,BinSize,0,max(r(1,:)));%First raster for state 1
    raster{2,1} = [];%Second raster for state 1
    raster{3,1} = [];%Third raster for state 1

    % Compute the Correlation matrices and the mean from the rasters, 
    datalenRaw = zeros(3,1);
    for i=1:size(raster,1)
        mAv{i,BinIndex} = zeros(N,1);
        CAv{i,BinIndex} = zeros(N,N);
        C1Av{i,BinIndex} = zeros(N,N);
        for j=1:size(raster,2)
            if ~isempty(raster{i,j})
                datalenRaw(i,j) = size(raster{i,j},2);
                mRaw{i,j} = mean(raster{i,j},2);
                CRaw{i,j} = raster{i,j}*transpose(raster{i,j})/datalenRaw(i,j);
                C1Raw{i,j} = raster{i,j}(:,1:(datalenRaw(i,j)-1))*transpose(raster{i,j}(:,2:datalenRaw(i,j)))/(datalenRaw(i,j)-1);                            
                mAv{i,BinIndex} = mAv{i,BinIndex} + datalenRaw(i,j)*mRaw{i,j};
                CAv{i,BinIndex} = CAv{i,BinIndex} + datalenRaw(i,j)*CRaw{i,j};
                C1Av{i,BinIndex} = C1Av{i,BinIndex} + datalenRaw(i,j)*C1Raw{i,j};
            end
        end
        datalen(i) = sum(datalenRaw(i,:));
        mAv{i,BinIndex} = mAv{i,BinIndex}/datalen(i);
        CAv{i,BinIndex} = CAv{i,BinIndex}/datalen(i);
        C1Av{i,BinIndex} = C1Av{i,BinIndex}/datalen(i);
    end

    % Computing for each state the distribution Pth and Pexp
    for state=1:1%Here there is only one state, but you can have more by just modifying the loading procedure abov. 
        for j=1:size(raster,2)
            rast{j} = raster{state,j};
        end
        m=mAv{state,BinIndex};
        C=CAv{state,BinIndex};
        C1=C1Av{state,BinIndex};
        [mtot,Ctot] = ConvertSpatial(m,C,C1);
        % Tanaka first guess
        [htotguess,Jtotguess]=hJinitial2(mtot,Ctot,2*N);
        [hstatGuess,JstatGuess]=hJinitial2(m,C,N);       
        % Monte-Carlo gradient descent
        [htot,Jtot,time_diff{state,BinIndex}]=GradientMC(mtot,Ctot,htotguess,Jtotguess,1000000,0.1,200,0.005);
        [hstat,Jstat,time_diff_stat{state,BinIndex}]=GradientMC(m,C,hstatGuess,JstatGuess,1000000,0.1,100,0.005);  
        % Storing the corresponding coefficients for the different bin
        % sizes and states
        [h,J,J1,hr,Jr] = DeConvertSpatial(htot,Jtot,hstat,Jstat);
        coefsh(:,state,BinIndex) = h;
        coefsJ(:,:,state,BinIndex) = J;
        coefsJ1(:,:,state,BinIndex) = J1;
        coefshr(:,state,BinIndex) = hr;
        coefsJr(:,:,state,BinIndex) = Jr;
        coefshstat(:,state,BinIndex) = hstat;
        coefsJstat(:,:,state,BinIndex) = Jstat;
    
        % Loop on the template size
        for TempSize=2:SizeMax
            
            %%% Model with spatio-temporal correlations %%%
            [Pexp,Pth,Features,Plist,RasterSize]=OcPred(rast,hstat,Jstat,hr,Jr,h,J,J1,TempSize,nbpatterns);%Pick some patterns, estimate empirically and compare o the prediction of the model. 
            divjs(BinIndex,TempSize,state)=PlotOcFit(Pexp,Pth,Features,BinSize,TempSize,N,RasterSize);% Plot the result
            if (PrintFig==1)
                print('-depsc',[DataFigures,'Data_s',int2str(state), '_b',int2str(BinSize),'_' int2str(N) 'x',int2str(TempSize),'.eps']);
            end
            PexpTab(1:length(Pexp),BinIndex,TempSize,state) = Pexp;
            PthTab(1:length(Pexp),BinIndex,TempSize,state) = Pth;
            PlistTab{BinIndex,TempSize,state} = Plist;
            FeaturesTab(1:length(Pexp),BinIndex,TempSize,state) = Features;
            RasterSizeTab(state,BinIndex) = RasterSize;
            ['Divergence, s' int2str(state) ', b ' int2str(BinSize) ' t ' int2str(TempSize),': ', num2str(divjs(1,BinIndex,TempSize,state)) ]
            close;
            
            %%% Model with spatial correlations %%%
            PthNoJ1 = OcPredFast(Plist,Pexp,hstat,Jstat,zeros(N,1),zeros(N,N),hstat,Jstat,zeros(N,N));
            divjsNoJ1(BinIndex,TempSize,state)=PlotOcFit(Pexp,PthNoJ1,Features,BinSize,TempSize,N,RasterSize);
            if (PrintFig==1)
                print('-depsc',[DataFigures,'DataNoJ1_s',int2str(state), '_b',int2str(BinSize),'_' int2str(N) 'x',int2str(TempSize),'.eps']);
            end
            PthNoJ1Tab(1:length(Pexp),BinIndex,TempSize,state) = PthNoJ1;
            ['Divergence when no Temp, s' int2str(state) ', b ' int2str(BinSize) ' t ' int2str(TempSize),': ', num2str(divjsNoJ1(1,BinIndex,TempSize,state)) ]
            close;
        
            %%% Model without correlation %%%
            PthNoJ = OcPredFast(Plist,Pexp,hstat,zeros(N,N),zeros(N,1),zeros(N,N),hstat,zeros(N,N),zeros(N,N));
            divjsNoJ(BinIndex,TempSize,state)=PlotOcFit(Pexp,PthNoJ,Features,BinSize,TempSize,N,RasterSize);
            if (PrintFig==1)
                print('-depsc',[DataFigures,'DataNoJ_s',int2str(state), '_b',int2str(BinSize),'_' int2str(N) 'x',int2str(TempSize),'.eps']);
            end
            PthNoJTab(1:length(Pexp),BinIndex,TempSize,state) = PthNoJ;
            ['Divergence when no J, s' int2str(state) ', b ' int2str(BinSize) ' t ' int2str(TempSize),': ', num2str(divjsNoJ(1,BinIndex,TempSize,state)) ]
            close;

        end
    end
end

save([pwd,'/Workspace/DataIsing.mat'])
