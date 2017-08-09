function [DATA] = getRowsContainingAValue(DATA, column)
%Function [DATA] = getRowsContainingAValue(DATA, column)
%
% Returns matrix DATA with rows that contain some value in specified column. 

rowcount = rowCount(DATA);
colcount = columnCount(DATA);

disp(['Getting rows that contain a value in column ' num2str(column) '.']);
disp(['Datamatrix contains ' num2str(rowcount) ' rows before operation.']);

empty_rows = [];

datacolumn = DATA{column};

for i=1:rowcount
    % check if the value on the column was string or numerical and make the
    % comparison according to that
    testid = datacolumn(i);
    
    if iscell(testid)
        testid = testid{1};
    end
        
    if isempty(testid)
        empty_rows = [empty_rows i];
    end
end

% remove every found value
for j=1:colcount
   DATA{j}(empty_rows) = [];
end

rowcount = rowCount(DATA);
disp(['Datamatrix contains ' num2str(rowcount) ' rows after operation.']);
