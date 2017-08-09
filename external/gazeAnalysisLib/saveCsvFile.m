function saveCsvFile(filename, headers, varargin)
%Function saveCsvFile(filename, headers, varargin)
%
% Function saves the data to a csv-file. Varargin must contain as many
% columns as headers contains column headers. A new function is made
% because Matlab writecsv contains some problematic aspects, e.g. in
% systems with different separator-symbol.

disp(['Writing the results to csv-file ' filename '...']);

fid = fopen(filename, 'w');

if fid == -1
   disp('Could not open the file for writing.');
   return;
end

% print first 2 lines
fprintf(fid, 'sep=,\n');
fprintf(fid, '%s', headers{1});
for i=2:length(headers)
    fprintf(fid, ',%s', headers{i});
end
fprintf(fid,'\n');

rowcount = length(varargin{1});
colcount = length(headers);
for i=1:rowcount
    if iscell(varargin{1})
        fprintf(fid, '%s', varargin{1}{i});
    else
        fprintf(fid, '%s', num2str(varargin{1}(i)));
    end
    
    for j=2:colcount
        if iscell(varargin{j})
            fprintf(fid, ',%s', varargin{j}{i});
        else
            fprintf(fid, ',%s', num2str(varargin{j}(i)));
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);

disp('Done.');