function [isgood] = testDataConsistency(DATA)
%Function [isgood] = testDataConsistency(DATA)
%
% This function performs tests to see if the data is in right form (e.g.
% all the columns have same amount of rows). Additional tests may be
% included later on.

rowcount = rowCount(DATA);
colcount = columnCount(DATA);

isgood = 1;

for i=1:colcount
    % if all the columns have same amount of rows
    if length(DATA{i}) ~= rowcount
        isgood = 0;
    end
end

if isgood
   disp('Data consistency test passed.');
else
   disp('Data consistency test not passed.');
end