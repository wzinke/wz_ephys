%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

function Init_dataset(name, pots_in, nx, ny, dx, dy, dt, t0, comment)

datasetdir = 'datasets';

VX = 1:0.05:nx;
VY = 1:0.05*(dx/dy):ny;

potsname = fullfile(datasetdir, [name '-pots.mat']); %#ok<NASGU>

if ~exist(potsname, 'file')
    pots = reshape(pots_in, [], nx, ny);
    missing = squeeze( abs(mean(pots,1)) + abs(std(pots,1)) )==0;
    % Zero mean and standard dev. means a missing signal
    n_miss = sum(missing(:));
    if n_miss>0
        disp(['  ' num2str(n_miss) ' missing signals detected!'])
        disp('  Patching data with local averages.')
        [msx, msy] = find(missing==1);
        for ii=1:n_miss
            ngs = neighbours([nx ny], msx(ii), msy(ii));
            pots(:,msx(ii),msy(ii)) = mean(pots_in(:,ngs),2);
        end;
    else
        msx = []; %#ok<NASGU>
        msy = []; %#ok<NASGU>
    end;
    disp('  Interpolating voltages')
    pots2 = interp2d(pots, VX, VY, 'splinem'); %#ok<NASGU>
    save(potsname, 'pots', 'pots2', 'name', 'dt', 't0', ...
        'dx', 'dy', 'nx', 'ny', 'msx', 'msy', 'comment');
    Get_csd(name, 'trad', 'no');
else
    disp('Dataset already exists. ')
    % What to do now?
end;
disp('  Done.')

function ngs = neighbours(sz,ii,jj) % 2 dim

ngs = [];
nnx = sz(1);
nny = sz(2);

if ii<nnx
    nr = sub2ind(sz, ii+1, jj);
    ngs = [ngs nr];
end;

if jj<nny
    nr = sub2ind(sz, ii, jj+1);
    ngs = [ngs nr];
end;

if ii>1
    nr = sub2ind(sz, ii-1, jj);
    ngs = [ngs nr];
end;

if jj>1
    nr = sub2ind(sz, ii, jj-1);
    ngs = [ngs nr];
end;
