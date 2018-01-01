function add_all(cwd)
% adds the specified directory to the path including all subdirectories.

%% Get the directory where this script resides to add the paths relatively to it
[pathStr,~,~] = fileparts(mfilename('fullpath'));

% generate all paths from this root
a = genpath(pathStr);
b = textscan(a,'%s','delimiter',':');
b = b{1};

% remove git directories
b(~cellfun(@isempty,strfind(b,'.git'))) = [];

% don't add external directory
b(~cellfun(@isempty,strfind(b,'/external/'))) = [];

% add all paths
addpath(b{:})

display([pathStr, ' added to the path']);

%  
