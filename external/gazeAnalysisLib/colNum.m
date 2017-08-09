function [column_number] = colNum(HEADERS, column_id)
%Function [column_number] = colNum(HEADERS, column_id)
%
% Finds the order number of the column_id and returns it as an integer. If
% column of that name not found, returns -1.

column_number = -1;
found = 0;

for n=1:length(HEADERS)
    if strcmp(HEADERS{n}{1}, column_id) && ~found
        column_number = n;
    end
end