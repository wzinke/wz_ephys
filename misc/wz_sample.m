function ov = wz_sample(vec, ns, re, nr)
% take <ns> random samples from input vector, <vec> with replacement (<re> = 1)
% or without (<re> = 0, default). 
% If <ns> is equal or larger than length(vec)). If only <vec> is provided 
% a bootstrap sample will be returned, i.e. a vector of the same length as <vec>
% with samples drawn with replacement (<re> = 1).
%
% <nr> specifies the number of samples drawn from <vec>. Note. if ns has
% two elements, the second one will be used as <nr>.
%
% Currently, it is assumed that the input is a row vector. 
% Thus, if nr > 1 the rows are the the <nr> samples. 
%
% This function is inspired by the sample function of the R package.
% (ok, Matlab provides now the randsample function that has a similar functionality but much slower.)
%
% ToDo: use a weighted sampling by allowing another vector defining the probability of each element.
%
% wolf zinke, 12.1.2014

lv = length(vec);

%%% how many samples?
if(exist('ns','var') == 0 || isempty(ns) == 1)
    ns = lv;
end

%%% replacement or not?
if(exist('re','var') == 0 || isempty(re) == 1)
    if(ns >= lv)
        re = 1;
    else
        re = 0;
    end
elseif(re == 0 && ns > lv)
    warning('Can not produce more data without replacement! \n Output will be generated with replacement instead.');
    re = 1;    
end

% how many replications?
if(exist('nr','var') == 0 || isempty(nr) == 1)
    if(length(ns) ==2)
        nr = ns(2);
    else
        nr = 1;
    end
end

%%% take the samples now!
if(re == 1)
    %%% with replacement
    ov = vec(randi(lv, nr, ns(1)));
else
    %%% without replacement�
%    elmat = repmat(1:lv,nr,lv);
    % nasty solution, have not found an option to just shuffle a 2D matrix.
    ov = nan(nr,ns);   
    for(i=1:nr)
        ov(i,:) = vec(randperm(lv,ns));
    end   
end
