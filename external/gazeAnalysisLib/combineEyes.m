function [newdatacol, newvalcol] = combineEyes(DATA, col1, col2, c1val, c2val, accepted_validities)
%Function [newdatacol, newvalcol] = combineEyes(DATA, col1, col2, c1val, c2val, accepted_validities)
%
% combines two columns to one by using the validity of columns. If both
% columns have validity among those specified in parameter
% accepted_validities, the mean of both columns is used as the value in
% new column and the validity stamp is taken from column 1. If either of
% the eyes has bad validity for certain datapoint and the other eye is fine,
% the good eye is used and the validity from that eye. If both eyes are
% bad, -1 is used as eye value and validity value. colX, cXval-parameters 
% are column numbers. accepted_validities is a cell array of validities
% considered "good". newdatacol  and newvalcol are vectors formed by the
% function containing result-coordinate and result-validity.

disp('Combining two columns into one.');

badeyec = -1;
badval = -1;

newdatacol = zeros(rowCount(DATA),1);
newvalcol = zeros(rowCount(DATA),1);


% for each datapoint
for i=1:rowCount(DATA)
    
    
    if ismember(DATA{c1val}(i), accepted_validities) && ismember(DATA{c2val}(i), accepted_validities)
        % both validities are okay
        newdatacol(i) = mean([DATA{col1}(i), DATA{col2}(i)]);
        newvalcol(i) = DATA{c1val}(i);
        
    elseif ismember(DATA{c1val}(i), accepted_validities) && ~ismember(DATA{c2val}(i), accepted_validities)
        newdatacol(i) = DATA{col1}(i);
        newvalcol(i) = DATA{c1val}(i);
    
    elseif ismember(DATA{c2val}(i), accepted_validities) && ~ismember(DATA{c1val}(i), accepted_validities)
        newdatacol(i) = DATA{col2}(i);
        newvalcol(i) = DATA{c2val}(i);
    
    else
        newdatacol(i) = badeyec;
        newvalcol(i) = badval;
        
    end
end

disp('Done.');