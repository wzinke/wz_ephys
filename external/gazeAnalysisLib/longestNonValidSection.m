function [longest_streak] = longestNonValidSection(DATA, valcol, timecol, accepted_validities)
%Function [longest_streak] = longestNonValidSection(DATA, valcol, timecol, accepted_validities)
%
% Calculates the longest time when value in validity-column is not part of 
% accepted_validities. The result is calculated according to the 
% time in timecol. valcol, timecol are column numbers while
% accepted_validities is a cell-array of validity-markings considered
% acceptable.

disp('Calculating longest non-valid section...');

datavector = DATA{valcol};
timevector = DATA{timecol};

longest_streak = 0;

if length(timevector)<2
   return; 
end

current_start = timevector(1);

for i=2:length(datavector)
    
    % if current streaks continue?
    if ~ismember(datavector(i), accepted_validities) % invalid data
        current = timevector(i) - current_start;
    else
        current_start = timevector(i);
        current = 0;
    end
    
    % if current streak surpasses longest one?
    if current > longest_streak
        longest_streak = current;
    end
    
end


disp('Done.');
