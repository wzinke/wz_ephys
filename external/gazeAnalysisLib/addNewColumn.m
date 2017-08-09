function [DATA, HEADERS] = addNewColumn(DATA, HEADERS, newcolumn, newcolumntitle)
%Function [DATA, HEADERS] = addNewColumn(DATA, HEADERS, newcolumn, newcolumntitle):
%
% Function adds new column to the DATA matrix and HEADERS-cell vector.
% After this operation the column can be operated as it was always part of
% .gazedata.

DATA{end+1} = newcolumn;
HEADERS{end+1} = {newcolumntitle};