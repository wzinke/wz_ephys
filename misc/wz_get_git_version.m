function [gitversion, moddate, GiNFO] = wz_get_git_version(gitdir, flnm)
% get information about the git version

% BUG: Con not cope with the existence of two occurrences of the same filename at
% different subdirectories of the Git repository...

moddate = [];
GiNFO   = [];

if(~exist('flnm', 'var'))
    flnm = [];
else
    [~, flstem, flext] = fileparts(flnm);
    flnm  = [flstem, flext];
end

if(exist(gitdir) == 0)
    error('File or directory not found!')
elseif(exist(gitdir,'file') == 2)
    [gitdir, flstem, flext] = fileparts(gitdir);
    flnm  = [flstem, flext];
elseif(exist(gitdir,'dir'))
    gitdir = fileparts(mfilename('fullpath'));
end

if(isempty(gitdir))
    gitdir = pwd;
end

cwd = pwd;
cd(gitdir);

[~, gitroot] = system('git rev-parse --show-toplevel');
GiNFO.gitroot = gitroot(1:end-1);
cd(GiNFO.gitroot);

if(nargin > 1 && ~exist(fullfile(gitdir,flnm), 'file'))
    [err, gitfl] = system(['git log --pretty=format: --name-status | cut -f2- | sort -u | grep ', flnm]);
    if(isempty(gitfl) || err)
        error('File not found!')
    end
end

[~, ingit] = system('git rev-parse --git-dir 2> /dev/null;');

if(isempty(ingit))
    cd(cwd);
    error('No git repository found!');
end

%[got_it, gitversion] = system('git describe --all --long | cut -d "-" -f 3 ');

if(~isempty(flnm))
    [err, gitfl] = system(['git log --pretty=format: --name-status | cut -f2- | sort -u | grep ', flnm]);
    gitfl(end) = [];
    GiNFO.filename  = gitfl;

    [got_it, gitversion] = system(['git log --follow  --pretty=format:''%h'' -- ', gitfl, ' | head -1']);
    [~, moddate]         = system(['git log -1 --format=%cd -- ', gitfl, ' | cut -d- -f 1 | cut -d'' '' -f2-']);
    [~, GiNFO.hash]      = system(['git hash-object ', gitfl]);
    [~, GiNFO.revcnt]    = system(['git log --follow  --pretty=format:''%h'' -- ', gitfl, ' | wc -l']);

    [got_it, GiNFO.repoversion] = system('git rev-parse --short HEAD');
    [~, GiNFO.repohash]         = system('git rev-parse HEAD');
    [~, GiNFO.reporevcnt]       = system('git rev-list HEAD --count');

    GiNFO.repoversion = [];
    GiNFO.repohash    = [];
    GiNFO.reporevcnt  = [];
else
    [got_it, gitversion] = system('git rev-parse --short HEAD');
    [~, moddate]         = system('git log -1 | grep Date | cut -d: -f2- | cut -d- -f1 | cut -d'' '' -f5-');
    [~, GiNFO.hash]      = system('git rev-parse HEAD');
    [~, GiNFO.revcnt]    = system('git rev-list HEAD --count');
end

moddate(end-1:end) = [];
gitversion(end)    = [];
GiNFO.gitversion   = gitversion;
GiNFO.date         = moddate;
GiNFO.hash(end)    = [];
GiNFO.revcnt(end)  = [];
GiNFO.timestamp    = datenum(datevec(GiNFO.date,'mmm dd HH:MM:SS yyyy'));

if(got_it ~= 0)
    warning('No information about git version available');
    gitversion = [];
    moddate    = [];
end

if(~isempty(gitdir))
    cd(cwd);
end
