function [ISIcell] = InterSpikeInformation (CellsInterested, DataCell)
    
ISIcell = cell(size(DataCell));

for j = 1 : length(CellsInterested)
  
    for i = CellsInterested(j)
        
        %Inter Spike Intervals  
        ISIcell{i , 1} = diff(DataCell{i , 3});   
        
        if isempty(ISIcell{i , 1}) == 0;    
            %Maximum inter Spike Interval   
            ISIcell{i , 2} = max(ISIcell{i , 1});    
            %Minimum inter Spike Interval  
            ISIcell{i , 3} = min(ISIcell{i , 1});   
            %Mean inter Spike Interval   
            ISIcell{i , 4} = mean(ISIcell{i , 1});    
            %STD inter Spike Interval   
            ISIcell{i , 5} = std(ISIcell{i , 1});   
            %Median inter Spike Interval  
            ISIcell{i , 6} = median(ISIcell{i , 1});    
            %ISI Histogram from 0 to 20000 ms    
            ISIcell{i , 7} = histc(ISIcell{i , 1} , (0:1:20000));
            %Cumulative Histogram   
            ISIcell{i , 8} = cumsum(ISIcell{i , 7});  
            %CMA  
            ISIcell{i , 9} = ISIcell{i , 8} ./ (1:1:20001); 
            %SkewnessCMA
            if  ISIcell{i , 2} > 20000
            ISIcell{i, 10} =  skewness(ISIcell{i , 9});
            else
                ISIcell{i, 10} = skewness(ISIcell{i , 9}(1: ISIcell{i , 2}));
            end
            %Skewness

                ISIcell{i, 11} = skewness(ISIcell{i , 1});%(1: ISIcell{i , 2}));

            
        end
        
    end
    
end         