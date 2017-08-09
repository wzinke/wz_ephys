function [Pexp,Pth,Feature,PatternList,RasterSize]=OcPred(rast,hstat,Jstat,hr,Jr,h,J,J1,TempSize,nbpatterns)
%Pick up a certain number of patterns and estimate both the empirical and
%model-predicted estimation
%inputs:
%rast: the raster
%h and J parameters of the model
%TempSize: size of the spatio-temporal patterns to be picked-up. 
%nbpatterns: number of spatio-temporal patterns to be picked-up. 
%
%Outputs:
%Pexp: empirically-estimated probabilities of the patterns
%Pth: theoretical prediction of the probabilities
%PatternList: list of the patterns picked-up
%RasterSize: effective length over which the patterns have been searched. 

N=length(h);

%We concatenate the different rasters and keep an index table to eliminate
%false detections at the borders of the rasters
raster = rast{1};
len= size(rast{1},2);
fakeind = ((len-TempSize+2):len);

for j=2:size(rast)
    raster = [raster rast{j} ];
    len = len + size(rast{j},2);
    fakeind = [fakeind ((len-TempSize+2):len)];
end

lendata=size(raster,2);

if (nargin<6)
    nbpatterns=floor(lendata/TempSize);
end
RasterSize = len-length(fakeind);

PatternList=zeros(N,TempSize,nbpatterns);
ind=1;
step=1;%index of the current potential pattern
while ((ind<=nbpatterns) && ((step+TempSize-1)<=size(raster,2)))
    pattern=raster(:,step:(step+TempSize-1));%We take a new pattern
    test1=pattern'*raster/N;
    test2=zeros(1,lendata-TempSize+1);
    for k=1:TempSize
        test2=test2+test1(k,k:(lendata-TempSize+k));
    end
    
    if (length(find((test2(1:(step-1))/TempSize)==1))==0)&&(isempty(find(fakeind==step)))
        %The pattern was not here before and is not at the borders of the rasters
        test2(fakeind)=0;%To avoid false detections at the borders
        Pexp(ind) = length(find((test2/TempSize)==1));
        Pexp(ind)=Pexp(ind)/(lendata-TempSize+1);
        %Compute the associated feature
        Feature(ind)=FeaturePattern(pattern);
        %We also store the pattern in a table (can be disabled if too slow):
        if nargin>3
            PatternList(:,:,ind)=pattern(:,:);
        end
    	%Then compute the theoretical prediction
        Eth(ind)=hstat'*pattern(:,1)+0.5*pattern(:,1)'*Jstat*pattern(:,1);
        for i=2:TempSize
            Eth(ind)=Eth(ind)+dot(pattern(:,i-1),J1*pattern(:,i)) + h'*pattern(:,i)+0.5*pattern(:,i)'*J*pattern(:,i) +  hr'*pattern(:,i-1)+0.5*pattern(:,i-1)'*Jr*pattern(:,i-1);
        end
        Pth(ind)=exp(Eth(ind));
        ind=ind+1;
    end
    step=step+1;
end

if (ind<nbpatterns)
    ['Number of patterns: ' int2str(ind-1)]
    PatternListTemp = PatternList(:,:,1:(ind-1));
    PatternList = PatternListTemp;
end


Pth=Pth*sum(Pexp)/sum(Pth);%Normalisation, faster than estimating the Z, and makes not that much difference
