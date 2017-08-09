function c = CSI()

global clustattrib;
global clustdata;

burstlength = .01;
overlap = [];
CSIval = [];
for i = 1:length(clustattrib.clusters)
      index = [];
      try   
         index = clustattrib.clusters{i}.index;   
      end
      for j = i+1:length(clustattrib.clusters)
         if (j~=i)
            index2 = [];
            try
               index2 = clustattrib.clusters{j}.index;   
            end
            if ((~isempty(index))&(~isempty(index2)))
               hitcount = 0;
               for k = 1:length(index)
                  
                  
                  time1 = clustdata.params(index(k),1)/10000;
                  amps = clustdata.params(index(k),2:5);
                  [maxamp, maxchannel] = max(amps);
                  
                  [timemin, timeindex] = min(abs(time1-(clustdata.params(index2,1)/10000)));
                  
                     
                  time2 = clustdata.params(index2(timeindex),1)/10000;
                  amps2 = clustdata.params(index2(timeindex),2:5);
                  maxamp2 = amps2(maxchannel);
                  timediff = time1-time2;
                  ampdiff = maxamp-maxamp2;
                  if (((time1-time2)<0)&((time1-time2)>-burstlength)&(ampdiff>0))                   
                     hitcount = hitcount+1;
                  elseif (((time1-time2)>0)&((time1-time2)<burstlength)&(ampdiff<0))                   
                     hitcount = hitcount+1;
                  end
               end  
               CSIval(i,j) = hitcount/length(index);
                          
            end
         end
      end
end

CSIval

