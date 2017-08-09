function [DATA] = removeRowsContainingValue(DATA, column, value)
%Function [DATA] = removeRowsContainingValue(DATA, column, value)
%
% Takes the rows that do not contain the specified value in specified
% column. The rows containing the specified value are discarded.

rowcount = rowCount(DATA);
colcount = columnCount(DATA);

disp(['Datamatrix contains ' num2str(rowcount) ' rows before operation.']);

rows_found = [];

for i=1:rowcount
    % check if the value on the column was string or numerical and make the
    % comparison according to that
    testid = DATA{column}(i);
    if iscell(DATA{column}(i))
       cnts = strcmp(testid, value);
    else
       cnts = testid == value;
    end

    % if row contains that 
    if cnts
        rows_found = [rows_found i];
    end
end

% remove every bad value
for j=1:colcount
   DATA{j}(rows_found) = [];
end

rowcount = sampleCount(DATA);
disp(['Datamatrix contains ' num2str(rowcount) ' rows after operation.']);
