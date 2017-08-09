%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BATCH FOR GLAUBER MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

i2mPath;
SizeMax=5;
nbpatterns=1000;
BinSize = 10;
datalen=100000;
N=8;
state=1;
BinIndex = 1;
tauTab = [1.5 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initializing storage matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PexpTab = zeros(nbpatterns,SizeMax);
PthTab = zeros(nbpatterns,SizeMax);
PthNoJ1Tab = zeros(nbpatterns,SizeMax);
PthNoJTab = zeros(nbpatterns,SizeMax);
PthNoJ0Tab = zeros(nbpatterns,SizeMax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating the Glauber raster and statistics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%J1test = random('unif',-0.15,0.15,N,N);
J1test = -0.1 + 0.2*rand(N,N);%or: random('unif',-0.1,0.1,N,N);
htest = -1.05 + 0.05*rand(N,1);%or: random('unif',-1.05,-1.0,N,1);

for k=1:length(tauTab)%Here we generate raster with the Glauber model for two different taus
    %Alternatively, a raster can be loaded by using the LoadRaster and
    %BinRaster functions
    tau = tauTab(k);
    raster{k} = GlauberRaster(J1test,htest,N,datalen,1,tau);%For this parameter tau, we generate a raster
    
    %HERE START THE ANALYSIS WITH THE STATISTICAL MODEL. In this example, the data are in
    %raster{k}, which is a binned raster composed of -1 and 1. See the
    %sub-directory /LoadRaster to load from data instead of using the
    %Glauber model. 
    m{k} = mean(raster{k},2);
    C{k} = raster{k}*transpose(raster{k})/datalen;
    C1{k} = raster{k}(:,1:(datalen-1))*transpose(raster{k}(:,2:datalen))/(datalen-1);
    [mtot,Ctot] = ConvertSpatial(m{k},C{k},C1{k});
    % First guess with the Tanaka approximation
    [htotguess,Jtotguess]=hJinitial2(mtot,Ctot,2*N); 
    [hstatGuess,JstatGuess]=hJinitial2(m{k},C{k},N);       
    
    % Monte Carlo gradient descent
    [htot,Jtot]=GradientMC(mtot,Ctot,htotguess,Jtotguess,10000000,0.1,50,0.007); 
    [hstat,Jstat]=GradientMC(m{k},C{k},hstatGuess,JstatGuess,1000000,0.1,100,0.007);
    [h,J,J1,hr,Jr] = DeConvertSpatial(htot,Jtot,hstat,Jstat); 
    
    %Store the fitted coefficients h and J in some tables for further use. 
    coefsh(:,k,BinIndex) = h;
    coefsJ(:,:,k,BinIndex) = J;
    coefsJ1(:,:,k,BinIndex) = J1;
    coefshr(:,k,BinIndex) = hr;
    coefsJr(:,:,k,BinIndex) = Jr;
    coefshstat(:,k,BinIndex) = hstat;
	coefsJstat(:,:,k,BinIndex) = Jstat;
        
    %Estimating the 3 models performance: predicting probabilities of
    %individual patterns (figure 1)
    for TempSize=1:SizeMax%For patterns of size 1 to SizeMax
    
        %%% Model with spatio-temporal correlations %%%
        rast{1} = raster{k};
        [Pexp,Pth,Features,Plist,RasterSize]=OcPred(rast,hstat,Jstat,hr,Jr,h,J,J1,TempSize,nbpatterns);%See OcPred picks up some patterns, compute their ocurrence empirically, and predict it from the model. 
    	divjs(TempSize,k)=PlotOcFit(Pexp,Pth,Features,BinSize,TempSize,N);% Plot the prediction vs empirical estimations
        
        PexpTab(1:length(Pexp),BinIndex,TempSize,k) = Pexp;% These 3 lines store the results in table for subsequent use in Fig 1.
        PthTab(1:length(Pexp),BinIndex,TempSize,k) = Pth;
        PlistTab{BinIndex,TempSize,k} = Plist;
        FeaturesTab(1:length(Pexp),BinIndex,TempSize,k) = Features;
        ['Divergence, t ' int2str(TempSize),': ', num2str(divjs(1,TempSize,k)) ]
    	close;

        %%% Model with only spatial correlations %%%    
        PthNoJ1 = OcPredFast(Plist,Pexp,hstat,Jstat,zeros(N,1),zeros(N,N),hstat,Jstat,zeros(N,N));
    	divjsNoJ1(TempSize,k)=PlotOcFit(Pexp,PthNoJ1,Features,BinSize,TempSize,N);
    	PthNoJ1Tab(1:length(Pexp),BinIndex,TempSize,k) = PthNoJ1;
        ['Divergence when no Temp, t ' int2str(TempSize),': ', num2str(divjsNoJ1(1,TempSize,k)) ]
        close;

        %%% Model without correlation %%%
        PthNoJ = OcPredFast(Plist,Pexp,hstat,zeros(N,N),zeros(N,1),zeros(N,N),hstat,zeros(N,N),zeros(N,N));
        divjsNoJ(TempSize,k)=PlotOcFit(Pexp,PthNoJ,Features,BinSize,TempSize,N);
        PthNoJTab(1:length(Pexp),BinIndex,TempSize,k) = PthNoJ;
        ['Divergence when no J, t ' int2str(TempSize),': ', num2str(divjsNoJ(1,TempSize,k)) ]
        close;
    
    end

    
end


save([pwd,'/Workspace/GlauberIsingData.mat'])
