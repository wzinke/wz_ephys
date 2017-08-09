function [percentage] = validGazePercentage(DATA, column, accepted_validities)
%Function [percentage] = validGazePercentage(DATA, column, accepted_validities)
%
% Function examines data and returns the amount of valid data. Accepted
% validities contains the validity marks to be considered "good".

disp(['Calculating validity percentage for column ' num2str(column) '...']);

% calculate the count of good data (value in the column 0)
total = length(DATA{column});
good = length(DATA{column}(ismember(DATA{column}, accepted_validities)));

percentage = good/total;

disp('Done.');
