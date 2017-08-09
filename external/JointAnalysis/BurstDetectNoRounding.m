function [TypeOfSpike,BurstStarts,BurstEnds,BurstDurations,NumSpikesInBursts,numberBursts, spikes] = BurstDetectNoRounding(DataCell,  CellsInterested, BurstThreshold,TailThreshold)
spikes = (DataCell{CellsInterested , 3});%spike times;
TypeOfSpike = zeros(1,length(spikes));


if spikes(2) - spikes(1) < BurstThreshold
    TypeOfSpike(1) = 1; % 1 for bursts, 2 for tails , 0 for individuals
end


for j = 2 : length(spikes)-1
    
    if spikes(j) - spikes(j-1) < BurstThreshold | spikes(j+1) - spikes(j)< BurstThreshold
        
        TypeOfSpike(j) = 1;
        
    end
end

if spikes(length(spikes)) - spikes(length(spikes)-1) < BurstThreshold
    
    TypeOfSpike(length(spikes)) = 1;
    
end

%% omitting bursts having fewer spikes then desired (in this case 2 spikes)

%+++comment out if two spike bursts are wanted
minNOspikes = 3; % specify the min number of spikes in a burst to be
BreakPointAfter=[];
BreakPointBefore=[];
burstspikes = find(TypeOfSpike==1);
    for j = 1:length(burstspikes)-1
        if spikes(burstspikes(j+1))- spikes(burstspikes(j)) > BurstThreshold
           BreakPointAfter = [BreakPointAfter burstspikes(j)];
        end
    end
BreakPointAfter = [BreakPointAfter burstspikes(end)];

for j = 2:length(burstspikes)
        if spikes(burstspikes(j))- spikes(burstspikes(j-1)) > BurstThreshold
           BreakPointBefore = [BreakPointBefore burstspikes(j)];
        end
    end
BreakPointBefore = [burstspikes(1) BreakPointBefore];

ZeroIndexes = find(BreakPointAfter-BreakPointBefore+ 1<minNOspikes); % spike indexes to be nullified
for l=1:length(ZeroIndexes)
TypeOfSpike(BreakPointBefore(ZeroIndexes(l)):BreakPointAfter(ZeroIndexes(l))) = 0;
end
%+++make comment out till here++++

%% TAIL CALCULATION

if spikes(2) - spikes(1) <= TailThreshold & ...
    TypeOfSpike(1) ~= 1
    TypeOfSpike(1) = 2; % 1 for bursts, 2 for tails , 0 for individuals
end


for j = 2 : length(spikes)-1
    
    if ((spikes(j) - spikes(j-1) <= TailThreshold)...
            | (spikes(j+1) - spikes(j)<=TailThreshold))...
            & TypeOfSpike(j) ~= 1
        
        TypeOfSpike(j) = 2;
        
    end
end

if (spikes(length(spikes)) - spikes(length(spikes)-1) <= TailThreshold) ...
    & TypeOfSpike(length(spikes)) ~= 1
    TypeOfSpike(length(spikes)) = 2;
end


%% omitting tails away from bursts and merging rest with the bursts

tailspikes = find(TypeOfSpike==2);
controlForLoop = TypeOfSpike; 
resultOfmerging=zeros(1,length(TypeOfSpike));
%/loop till all the burst related tails merged and counted as bursts
while isequal(controlForLoop,resultOfmerging)==0
    controlForLoop = TypeOfSpike;
    for m = 1 : length(tailspikes)
        if (tailspikes(m)== 1)...
                & (spikes(tailspikes(m) + 1) - spikes(tailspikes(m)) <= TailThreshold...
                & TypeOfSpike(tailspikes(m) + 1)==1)
            TypeOfSpike(tailspikes(m)) = 1;
        elseif tailspikes(m)== length(TypeOfSpike)...
                &(spikes(tailspikes(m))- spikes(tailspikes(m)-1) <= TailThreshold...
                & TypeOfSpike(tailspikes(m) - 1)==1)
            TypeOfSpike(tailspikes(m)) = 1;           
        elseif (tailspikes(m)~= 1)& tailspikes(m)~= length(TypeOfSpike)
            if(spikes(tailspikes(m) + 1) - spikes(tailspikes(m)) <= TailThreshold...
                & TypeOfSpike(tailspikes(m) + 1)==1) | (spikes(tailspikes(m))...
                - spikes(tailspikes(m)-1) <= TailThreshold...
                & TypeOfSpike(tailspikes(m) - 1)==1)           
            TypeOfSpike(tailspikes(m)) = 1;
            end
        end
    end
resultOfmerging = TypeOfSpike;
tailspikes = find(TypeOfSpike==2);
end

nonbursts= find(TypeOfSpike ~= 1);
TypeOfSpike(nonbursts)=0;
        
%%Burst Starts  and Burst Ends
burstspikes = find(TypeOfSpike==1);
if isempty(burstspikes)~= 1
BreakPoint=[];
    for j = 1:length(burstspikes)-1
        if spikes(burstspikes(j+1))- spikes(burstspikes(j)) > TailThreshold
           BreakPoint = [BreakPoint burstspikes(j)];
        end
    end
BreakPointEnds = [BreakPoint burstspikes(end)];
BurstEnds = spikes(BreakPointEnds);

BreakPoint=[];
    for j = 2:length(burstspikes)
        if spikes(burstspikes(j))- spikes(burstspikes(j-1)) > TailThreshold
           BreakPoint = [BreakPoint burstspikes(j)];
        end
    end
BreakPointStarts = [burstspikes(1) BreakPoint];
BurstStarts = spikes(BreakPointStarts);
NumSpikesInBursts = BreakPointEnds - BreakPointStarts +1;
BurstDurations = BurstEnds - BurstStarts;
numberBursts = length(BurstStarts);   
else
   BurstStarts = '-';
   BurstEnds= '-';
   BurstDurations= '-';
   NumSpikesInBursts= '-';
   numberBursts= '-'; 
end
