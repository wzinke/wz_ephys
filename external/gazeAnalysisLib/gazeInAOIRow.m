function [rownumber] = gazeInAOIRow(data, xcol, ycol, aoicoord, type)
%Function [rownumber] = gazeInAOIRow(data, xcol, ycol, aoicoord, type)
%
% Returns the rownumber when gaze is (according to last/first) inside the area aoi. 
% If the user does not enter aoi, then -1 is returned. If user does not leave 
% the aoi, -1 is returned.
% GAZE = mean of two eyes, all times
% aoicoord = [xmin xmax ymin ymax] specified in the x and y pixel
% coordinates
% type = 'first' returns the first row when gaze is in the aoi
%        'last' returns the last row gaze was in the aoi
%        'firstleave' returns the first row of gaze leaving the aoi

disp(['Finding row gaze passes aoi border [' num2str(aoicoord(1)) ' ' num2str(aoicoord(2)) ' ' ...
      num2str(aoicoord(3)) ' ' num2str(aoicoord(4)) '] type:' type '...']);

% find the time when gaze goes out of the aoi
timefound = 0;

x = data{xcol};
y = data{ycol};

switch type
    
    case 'first'
        
        % if not found, return -1
        rownumber = -1;
        
        for i=1:length(x)
    
            % both coordinates inside the specific aoi
            if (aoicoord(1) < x(i)  && x(i) < aoicoord(2)) && (aoicoord(3) < y(i) && y(i) < aoicoord(4))
                if ~timefound
                    timefound = 1;
                    rownumber = i;
                end
            
            end
        end
        
    case 'last'
        
        % if not found, return -1
        rownumber = -1;
        
        for i=1:length(x)
    
            % both coordinates inside the specific aoi
            if (aoicoord(1) < x(i)  && x(i) < aoicoord(2)) && (aoicoord(3) < y(i) && y(i) < aoicoord(4))
                rownumber = i;             
            end
        end
        
    case 'firstleave'
        
        % if not found to leave aoi, return max number of rows
        rownumber = rowCount(data);
        
        for i=1:length(x)
    
            % coordinates outside the specific aoi
            if (aoicoord(1) > x(i) || x(i) > aoicoord(2) || aoicoord(3) > y(i) || y(i) > aoicoord(4)) && ~timefound
                rownumber = i;
                timefound = 1;
            end
        end
end

disp('Done.');


% for i=1:length(x)
%     
%     % both coordinates inside the specific aoi
%     if (aoicoord(1) < x(i)  && x(i) < aoicoord(2)) && (aoicoord(3) < y(i) && y(i) < aoicoord(4))
%         if strcmp(type, 'first')
%             if ~timefound
%                 timefound = 1;
%                 rownumber = i;
%             end
%         elseif strcmp(type, 'last')
%             rownumber = i;             
%         end
%     else
%         if strcmp(type, 'firstleave') && ~timefound
%             rownumber = i;
%             timefound = 1;
%         end
%     end
% end
