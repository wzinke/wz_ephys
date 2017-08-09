function [files] = findGazeFilesInFolder(folder, ending)
%Function [files] = findGazeFilesInFolder(folder)
%
% Returns the full pathnames of gazefiles in the folder. The gazefiles 
% ending is specified as a parameter. (e.g. '.gazedata' or '.csv')

disp(['Retrieving the gazefiles in the folder ' folder '...']);

all_files = dir(folder);

filecounter = 0;
files = {};

for i = 1:length(all_files)
    [a, b, c] = fileparts(all_files(i).name);
    if strcmp(c, ending)
        filecounter = filecounter+1;
        files{filecounter} = [folder a b c];
    end
end

if filecounter == 0
    disp('No gazefiles in the folder.');
else
    disp('Done.');
end