function [tabPopn,x] = CorrStatePopn(raster,delta,BinSize,fig,norm,colors)
%Compute a table of the cross correlation from the raster
if (nargin<6)
    colors{1} = [1 0 0];
    colors{2} = [0 1 0];
    colors{3} = [0 0 1];
end


for i=1:size(raster,1)
    for j=1:size(raster,2)       
        if ~isempty(raster{i,j})
            rasterPopn{i,j} = sum(raster{i,j},1);
        else
            rasterPopn{i,j} = [];
        end
    end
end

[tab,x] = CorrStateTab(rasterPopn,delta,BinSize,fig,norm,colors);
size(tab)
tabPopn = tab;

