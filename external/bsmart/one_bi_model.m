function one_bi_model(dat,order,startp,window)
% ONE_BI_MODEL One window for bivariate model
%
% Usage:
%   one_bi_model(dat,order,startp,window);
%
% Input(s):
%   dat: input file in matlab format;
%   order: model order
%   startp: start position of the window
%   window: window length
%
% Output(s):
%   in the bsmartroot/Onewindow_Coefficient directory

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 09:40:21$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% parse directory
cdir = pwd;                 % find current directory
p = mfilename('fullpath');
fdir = fileparts(p);        % find function directory
cd(fdir);                   % change current dir to function dir

% processing
si=size(dat);
ntrails=si(3);
nchannels=si(2);
start=startp;
channel2=nchannels;
channel=2;
trail=ntrails;
points=window;
save channel channel -ascii;
save trail trail -ascii;
save points points -ascii;
save order order -ascii;

dat=dat(start:start+window-1,:,:);
PathName=fullfile(pwd,'Onewindow_Coefficient');
coefile = [PathName,'\*'];
delete(coefile);            % delete all coefficients files in the directory

for i=1:(channel2-1)
    for j=(i+1):channel2
        dat1=dat(:,i,:);
        dat2=dat(:,j,:);
        dat3=cat(2,dat1,dat2);
        writedat('dataset.bin',dat3);
        % TODO: don't return '0', use 'waitebar'
        if ispc
            eval(['unix ' '(''' 'opssfull ' 'dataset.bin ' ' A ' 'Ve ' 'AIC' ''')'])
        else
            eval(['unix ' '(''' './opssfull ' 'dataset.bin ' ' A ' 'Ve ' 'AIC' ''')']);
        end%if
        % save and move files
        ii=num2str(i);jj=num2str(j);
        FileName=['AR_C_' ii 'and' jj];
        fA=fullfile(PathName,FileName);
        % eval(['!MOVE ' 'A ' fA ]);
        movefile('A',fA);
        FileName2=['AR_N_' ii 'and' jj];
        fVe=fullfile(PathName,FileName2);
        % eval(['!MOVE ' 'Ve ' fVe ]);
        movefile('Ve',fVe);

        delete AIC
        delete dataset.bin
    end
end

delete channel
delete trail
delete points
delete order

helpdlg('All model coefficient save in Onewindow_Coefficient directory',...
    'Completed');

% restore directory
cd(cdir);

end%one_bi_model

% [EOF]