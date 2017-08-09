function [DATA] = clipLastRows(DATA, cutrow)
%Function [DATA] = clipLastRows(DATA, cutrow)
%
% Returns the last numrows rows of DATA-matrix in same format.

rowcount = rowCount(DATA);
colcount = columnCount(DATA);

disp(['Picking rows after row number ' num2str(cutrow) '...'])
disp(['Datamatrix contains ' num2str(rowcount) ' rows before operation (1.st column).']);

% if there are more rows than minimum
if cutrow < rowcount
    % put all the columns after numrows as blank
    for i=1:colcount
        DATA{i}(1:cutrow) = [];
    end
else
    disp('The data contains less datapoints than the cutrow-argument. Not cutting anything.');
end

rowcount = rowCount(DATA);

disp(['Datamatrix contains ' num2str(rowcount) ' rows after operation.']);