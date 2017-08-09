function [tab,x] = CorrStateTab(raster,delta,BinSize,fig,norm,colors)
%Correlation table for different states and raster block. 
%Same arguments than for CorrTab, except than:
%-the raster is here a cell array
%-colors is a cell array which indicates, for the plotting option, the colors associated to the
%different states

if (nargin<6)
    colors{1} = [1 0 0];
    colors{2} = [0 1 0];
    colors{3} = [0 0 1];
end

nbneur = size(raster{1,1},1);

tab = zeros(size(raster,1),2*delta+1,nbneur,nbneur);

for i=1:size(raster,1)
    tabTemp2 = zeros(2*delta+1,nbneur,nbneur);
    for j=1:size(raster,2)       
        if ~isempty(raster{i,j})
            lendata(i,j) = size(raster{i,j},2);
            [tabTemp,x] = CorrTab(raster{i,j},delta,BinSize,0,norm);
            %tab(i,:,:,:) = tab(i,:,:,:) + lendata(i,j)*tabTemp(:,:,:);
            tabTemp2 = tabTemp2 + lendata(i,j)*tabTemp;
        end
    end    
    tab(i,:,:,:) = tabTemp2(:,:,:)/sum(lendata(i,:));%tab(i,:,:,:)/sum(lendata(i,:));
end

x=BinSize*(-delta:delta);

if (fig>0)
    lowB = min(min(min(min((tab)))));
    upB = max(max(max(max(tab))));
    figure;
    for i=1:nbneur
        for j=1:nbneur
            subplot(nbneur,nbneur,nbneur*(i-1)+j);
            for s=1:size(raster,1)
                plot(x,tab(s,:,i,j),'color',colors{s});
                if (fig==2)
                    set(gca,'YLim',[lowB upB]);
                else
                    set(gca,'YLim',[min(min(tab(:,:,i,j))) max(max(tab(:,:,i,j)))]);
                end
                hold on
            end
        end
    end
end
