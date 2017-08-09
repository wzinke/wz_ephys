function pmkdir(thisdir)
% create a directory. Check first, if it does exist, if not, check if
% parent directory does exist and create it if not.

if(thisdir(1) ~= '/')
    thisdir = fullfile(pwd,thisdir);
end

if(~exist(thisdir,'dir'))
    dirstem = fileparts(thisdir);
    if(~exist(dirstem,'dir'))
        pmkdir(dirstem);
    else
        mkdir(thisdir);
    end
end
