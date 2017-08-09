function bbz_save_function(varargin)
% BBZ_SAVE_FUNCTION Saves active m-file with time stamp from command line
%  
% DESCRIPTION 
% This function saves
%  
% SYNTAX 
% BBZ_SAVE_FUNCTION;
% BBZ_SAVE_FUNCTION('saveas');
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 08 Jul 2013 20:34:29 CDT by bramzandbelt 
% $Modified: Tue 09 Jul 2013 09:13:12 CDT by bram

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. VARIABLE HANDLING
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get currently active m-file and determine if this is a new file
% =========================================================================
docObj          = matlab.desktop.editor.getActive;
isNew           = strncmpi(docObj.Filename,'untitled',8);

% Identify user
% =========================================================================
if ismac || isunix
    user        = getenv('USER');
elseif ispc
    user        = getenv('name');
end

% Define a time string containing date, time, and timezone
% =========================================================================
if isempty(which('now.py'))
    error('now.py not found in path');
else
    [~,timeStr] = unix(['python ',which('now.py')]);
end

% Determine what to do
% =========================================================================
if nargin == 0 && isNew
    toDo        = 'savenew';
elseif nargin == 0 && ~isNew
    toDo        = 'save';
elseif nargin == 1
    toDo        = 'saveas';
else
    error('Incorrect number of input arguments.');
end

% #.# Specify file name
% =========================================================================
switch toDo
    case {'savenew','saveas'}
        [f,p]   = uiputfile('*.m');
        fName   = fullfile(p,f);
    case 'save'
        [p,f,e] = fileparts(docObj.Filename);
        fName   = fullfile(p,[f,e]);
end

% #.# Update line containing $Modified tag
% =========================================================================

% Convert char array to cell array
textCell        = matlab.desktop.editor.textToLines(docObj.Text);

% Identify row with $Modified tag
iRow            = cellfun(@(a) strncmpi(a,'% $Modified',11),textCell);

% Update tage
textCell{iRow}  = sprintf('%% $Modified: %s by %s',deblank(timeStr),user);

% Convert cell array back to char array
docObj.Text     =  matlab.desktop.editor.linesToText(textCell);

% Save MATLAB function
docObj.saveAs(fName);
