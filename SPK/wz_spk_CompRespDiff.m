function DifComp = wz_spk_CompRespDiff(unit1, unit2, timwin, mincons, pthr, nBoot, smplsz, xvec)
% wz_spk_CompRespDiff - compare time of response difference between two neurons.
%                       Apply a bootstrrap approach to check for co-variation.
%  
% DESCRIPTION 
%
%  
% SYNTAX 
%   TST = wz_spk_CompRespDiff()
%
%   Input:
%       
%   
%   Output:
%       
%
% REFERENCES 
%
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 06-Oct-2014 by wolf zinke
% $Modified: 

%% get 


% number of bootstrap repetitions
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1;
end

% number of bootstrap repetitions
if(~exist('smplsz','var') || isempty(smplsz))
    if(nBoot == 1)
        smplsz = [];
    else
        smplsz = 25;
    end
end

if(~exist('timwin','var') || isempty(timwin))
    timwin = [30, 300];
end

if(~exist('spontwin','var') || isempty(spontwin))
    spontwin = [-500, 0];
end

if(~exist('krnl','var') || isempty(krnl))
    krnl = 'gauss';
end

if(~exist('kw','var') || isempty(kw))
    kw = 10;
end


