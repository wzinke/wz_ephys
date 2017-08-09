function times = maketimesfile(matclustfile)

load(matclustfile)

timelist = clustdata.timefilternames;
times = [];
count = 0;
for i = 2:length(clustdata.timefilternames)
    dash = strfind(timelist{i},'-');
    if ~isempty(dash)
        count = count+1;
        spaces = strfind(timelist{i},' ');
        lastspace = spaces(end);
        firstspace = spaces(1);
        name = timelist{i}(firstspace:lastspace);
        name(strfind(name,' ')) = '';
        name = lower(name);
        starttime = timelist{i}(lastspace+1:dash-1);
        endtime = timelist{i}(dash+1:end);
        times{count} = [name,' ',starttime,' ',endtime];
    end
end

