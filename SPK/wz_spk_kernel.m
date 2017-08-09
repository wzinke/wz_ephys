function krnl = wz_spk_kernel(ktype, kwdth)
% Determine a kernel of type <ktype> with the width <kwdth> to be convolved
% with a spike train.
%
% ToDo: if needed add gamma and ex-gauss kernels.
%       use option to adjust position according to center of mass or peak location.
%
%
% wolf zinke, 13.1.2014

if(exist('ktype','var') == 0 || isempty(ktype) == 1)
    ktype = 'gauss';
end

if(exist('kwdth','var') == 0 || isempty(kwdth) == 1)
    kwdth = 10;
end

% Calculate values of the kernel
switch lower(ktype)
    case {'gauss', 'norm'}
        Half_BW = round(4*kwdth); % This is a bit ad hoc...
        x       = [-Half_BW : Half_BW];
        krnl    = (1/(sqrt(2*pi)*kwdth)) * exp(-1*((x.^2) / (2*kwdth^2)));
%        krnl = exp(-(x.^2)./2./kwdth^2)./kwdth./sqrt(2*pi);

    case 'tri'
        Half_BW = round(sqrt(6)*kwdth);
        x       = [-Half_BW : Half_BW];
        krnl    = (1/(6*kwdth))*(sqrt(6)*kwdth - abs(x));

    case {'box', 'rect'}
        Half_BW = round(sqrt(3)*kwdth);
        x       = -Half_BW : Half_BW;
        krnl    = ones(1,length(x));

    case {'psp', 'exp'}  % 'exp' is a bit misplaced here but has its reason from the function history...
        % Be careful, asymmetric kernel. Use this only if you want to
        % mimic the post-synaptic effect of the neuronal response. I
        % discourage the use of asymmetric kernels for a statistical
        % description of the response characteristics.
        if(length(kwdth) == 2)
            tg=kwdth(2);
            kwdth(2) = [];
        else
            tg = 1;
        end
        % Use kwdth = 20 and tg = 1 to be in line with:
        % K. G. Thompson, D. P. Hanes, N. P. Bichot, and J. D. Schall
        % Perceptual and motor processing stages identified in the activity of
        % macaque frontal eye field neurons during visual search
        % J Neurophysiol December 1, 1996 76:(6) 4040-4055

        Half_BW = ceil(kwdth*8);
        x = [0 : Half_BW];
        krnl = [ zeros(1,Half_BW) , (1-(exp(-(x./tg)))).*(exp(-(x./kwdth)))];

    case {'lognorm', 'ln'}
        % same remarks here as just before about asymmetric kernels...
        Half_BW = ceil(kwdth*8);
        x = [0 : Half_BW];

        % krnl = [zeros(1,Half_BW), exp(-0.5 * (log(x) ./ kwdth).^2) ./ (x .* sqrt(2*pi) .* kwdth)];
        krnl = [zeros(1,Half_BW), lognpdf(x,0,kwdth)];

    case 'expo'  % should be 'exp' but see 'psp' for reasons.
        % same remarks here as just before about asymmetric kernels...
        Half_BW = ceil(kwdth*8);
        x = [0 : Half_BW];
        krnl = [zeros(1,Half_BW), exppdf(x,kwdth)];

    otherwise
        error(['Kernel type "', ktype, '" is not implemented!']);
end

krnl = krnl - krnl(end);  % avoid a general step offset introduced by the kernel.
krnl(krnl<0) = 0;
krnl = krnl/sum(krnl);


