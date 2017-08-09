function [DATA] = clipMilliSecondsAfter(DATA, HEADERS, millisec)
%Function [DATA] = clipMilliSecondsAfter(DATA, HEADERS, millisec)
%
% Returns the millisec milliseconds after specified time of DATA-matrix in 
% same format.

rowcount = rowCount(DATA);
colcount = columnCount(DATA);

disp(['Picking data after ' num2str(millisec) ' milliseconds from data using TETTime...']);
disp(['Datamatrix contains ' num2str(rowcount) ' rows before operation.']);

TETTime = DATA{colNum(HEADERS, 'TETTime')};

millisec_at_start = TETTime(1);

% if there are more rows than the ms limit
if millisec_at_start + millisec < TETTime(rowcount)
    
    cutrow = find(TETTime > millisec_at_start + millisec, 1, 'first');
    
    % put all the columns after numrows as blank
    for i=1:colcount
        DATA{i}(1:cutrow) = [];
    end
else
    % empty all data
    for i=1:colcount
        DATA{i}(1:end) = [];
    end
    
    disp('The data is shorter than millisecond-limit. Returning empty gazedata matrix.');
end

rowcount = rowCount(DATA);

disp(['Datamatrix contains ' num2str(rowcount) ' rows after operation.']);


%TODO:: TEST THIS!!!!
%!!!!!!!!