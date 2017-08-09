function [DATA, HEADERS] = loadGazedataFile(file, columncount, dataformat)
%Function [DATA, HEADERS] = loadGazedataFile(file, columncount, dataformat)
%
% Loads the gaze file and returns all datalines and columns in DATA-array and 
% headers in HEADERS-cell-vector. Load function is sensitive in content
% format and you need to specify that as a parameter. The format follows
% Matlab's common format e.g. in textscan function (e.g. '%s %f %f %s' for
% data where first column contains strings, second and third float-numbers 
% and fourth again strings). Column headers are assumed to be in the same 
% file and similar numbers as columns in the actual data.

disp(['Reading file ' file '...']);

str = fileread(file);

disp('Done.');
disp('Replacing uncompatible marks on the file...');

% corrections to the unwanted marks on the file (some usual found in
% gazedata files)
str = strrep(str, '-1.#INF', '-1');
str = strrep(str, '-1.#IND', '-1');
str = strrep(str, '1.#INF', '-1');
str = strrep(str, '1.#IND', '-1');
str = strrep(str, '1.#QNAN', '-1');
str = strrep(str, '-1.#QNAN', '-1');

disp('Done.');
disp('Loading file to datastructs...');

% form header format string
headerformat = '%s';
for i=1:columncount-1
    headerformat = [headerformat ' %s'];
end

%Read header lines
HEADERS = textscan(str, headerformat, 1);

%Read data columns
DATA = textscan(str, dataformat, 'HeaderLines', 1, 'Delimiter', '\t');

disp('File loaded.');