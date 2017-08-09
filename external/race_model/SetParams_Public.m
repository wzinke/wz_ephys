%% SetParams_Public.m
% =========================================================================
%                               set parameters
% =========================================================================

meanGO = param(1);
stdGO = param(2);
modelName = ' ';

if modelNo == 0  % Independent Model
    meanSTOP = param(3); 
    stdSTOP = param(4); 
    BGo = 0;
    BStop = 0;
    time_delay = ceil(param(5));
    modelName = 'Independent';
end

if modelNo == 1  % Interactive Model
    meanSTOP = param(3); 
    stdSTOP = param(4); 
    BGo = param(5);
    BStop = param(6);
    time_delay = ceil(param(7));
    modelName = 'Interactive';
end

if time_delay < 0
    time_delay = 0;
end

BGo = abs(BGo);
BStop = abs(BStop);

% leakage terms
kGo = 0;
kStop = 0;

TrGO = 35;  % set to Monkey C's value

% output values once
if nt == 1 & iteration == 0
    modelNo
    meanGO       
    stdGO
    meanSTOP
    stdSTOP
    BGo
    BStop
    TrGO
    time_delay
end