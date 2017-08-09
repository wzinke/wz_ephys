function [tab,x] = CorrTab(r,delta,BinSize,fig,norm)
%compute and plot auto and cross correlation from a raster. Not normalized
%Inputs:
%r: raster, with 1 and -1
%delta: number of different delays in the computed cross corr
%BinSize: the raster is already binned, This param allows to convert the
%x-axis into ms. 
% fig: 0=no plot, 1=plot with each of the graph having its own scale, 2=plot with same scale for all the graphs 
%
%Outputs:
%tab(t,cell1,cell2) contains the correlation between cell 1 and 2, with
%delay x(t). x vector (unit: ms) can serve for plotting purposes

len=size(r,2);

if (norm==1)
    m=mean(r,2);
    sigma=std(r,1,2);
    tab1=r-m*ones(1,len);%substract the mean
end

for i=1:2*delta+1
    d=i-delta-1;
    if (norm==0)
        tab(i,:,:) = (r(:,max(1,d+1):min(len,len+d))+1)*((r(:,max(1,1-d): min(len,len-d))+1)')/(4*(len-abs(d))*0.001*BinSize);
    else
        tab(i,:,:) = (tab1(:,max(1,d+1):min(len,len+d))*tab1(:,max(1,1-d): min(len,len-d))'/(len-abs(d)))./(sigma*sigma');
    end
end

x=BinSize*(-delta:delta);
% d<0: 1:len+d , 1-d:len
% d>0: d+1:len , 1:len-d
% max(1,d+1):min(len,len+d) , max(1,1-d): min(len,len-d)

if (fig>0)
    nbneur=size(r,1);
    lowB = min(min(min(tab)));
    upB = max(max(max(tab)));
    figure;
    for i=1:nbneur
        for j=1:nbneur
            subplot(nbneur,nbneur,nbneur*(i-1)+j);
            plot(x,tab(:,i,j));
            if (fig==2)
                set(gca,'YLim',[lowB upB]);
            end
        end
    end
end
