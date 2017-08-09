function [DATA] = clipFirstMilliSeconds(DATA, HEADERS, millisec)
%Function [DATA] = clipFirstMilliSeconds(DATA, HEADERS, millisec)
%
% Returns the first millisec milliseconds of DATA-matrix in same format.

rowcount = rowCount(DATA);
colcount = columnCount(DATA);

disp(['Picking first ' num2str(millisec) ' milliseconds from data using TETTime...']);
disp(['Datamatrix contains ' num2str(rowcount) ' rows before operation.']);

TETTime = DATA{colNum(HEADERS, 'TETTime')};

millisec_at_start = TETTime(1);

if millisec <= 0
    % if cut nothing
    
    % put all the columns after numrows as blank
    for i=1:colcount
        DATA{i}(:) = [];
    end
    
elseif millisec_at_start + millisec < TETTime(rowcount)
    % if there are more rows than the ms limit
    cutrow = find(TETTime < millisec_at_start + millisec, 1, 'last');
    
    % put all the columns after numrows as blank
    for i=1:colcount
        DATA{i}(cutrow+1:rowcount) = [];
    end
else
    disp('The data contains less datapoint than the millisecond-count. Not cutting anything.');
end

rowcount = rowCount(DATA);

disp(['Datamatrix contains ' num2str(rowcount) ' rows after operation.']);


%TODO:: TEST THIS!!!!