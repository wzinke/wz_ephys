function out = timetrans(times,UnitsPerSec,dir)

%out = timestrans(times, UnitsPerSec, dir)

%transforms time inputs from seconds to hours:minutes:seconds or vice versa

%TIMES is times that you want to transform, ie [60 3600] seconds or
%[0:01:00 1:00:00] in hours:min:secs
%UnitsPerSec is units in each sec
%dir is 1 if going from seconds to hours:min:sec, 2 if going from
%hrs:min:sec to sec

if (dir==1)
	
	
	
	for i = 1:length(times)
      if (times(i) > 0)
            hours = floor(times(i)/(60*60*UnitsPerSec));
	        minutes = floor(times(i)/(60*UnitsPerSec))-(hours*60);
	        seconds = floor(times(i)/(UnitsPerSec))-(hours*60*60)-(minutes*60);
      else
            hours = abs(ceil(times(i)/(60*60*UnitsPerSec)));
	        minutes = abs(ceil(times(i)/(60*UnitsPerSec))+(hours*60));
	        seconds = abs(ceil(times(i)/(UnitsPerSec))+(hours*60*60)+(minutes*60));
      end
      if (minutes<10)
          tempmin = ['0',num2str(minutes)];
      else
          tempmin = [num2str(minutes)];
      end
      if (seconds<10)
          tempseconds = ['0',num2str(seconds)];
      else
          tempseconds = [num2str(seconds)];
      end
      out{i,1} = [num2str(hours),':',tempmin,':',tempseconds];
      if (times(i)<0)
          out{i,1} = ['-',out{i,1}];
      end
    end
    
elseif (dir==2)
    for i = 1:length(times)
        t = [0 0 0 0 0 0];
        temptime = times{i};
        colons = findstr(temptime,':');
        startind = length(temptime);
        count = 6;
        for j = length(colons):-1:1
            t(count) = str2num(temptime((colons(j)+1):startind));
            startind = colons(j)-1;
            count = count-1;    
        end
        t(count) = str2num(temptime(1:colons(1)-1));
        out(i,1) = etime(t,[0 0 0 0 0 0])*UnitsPerSec;
    end
end
                
        