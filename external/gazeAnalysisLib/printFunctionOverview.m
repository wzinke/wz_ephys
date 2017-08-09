function printFunctionOverview
%Function printFunctionOverview
%
% Prints the overview of functions contained in the toolbox by displaying
% their help-texts.

% register the rootdirectory where the function is located
rootdir = fileparts(mfilename('fullpath'));

files = dir(rootdir);

for i=1:length(files)
    
    filename = files(i).name;
    
    [a, b, c] = fileparts(filename);
    
    if ~files(i).isdir && strcmp(c, '.m') && exist(filename) == 2
     %   disp(['Function ' filename ':']);

        help(filename)
        
        disp(' ');
        %disp(' ');
    end
end