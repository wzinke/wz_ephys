function add_all(cwd)
% adds the specified directory to the path including all subdirectories.

if(exist('cwd','var') == 0)
   cwd = pwd;
elseif(isempty(cwd) == 1)
   cwd = pwd;
end

addpath(genpath(cwd));


