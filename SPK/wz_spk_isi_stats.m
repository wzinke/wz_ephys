function [CV, mISI, SK] = SPK_isi_stats(isi_lst, DIM)
% Calculation of standard measurement of an ISI distribution such as the mean ISI,
% the coefficent of variation (CV) and the coefficent of the asymmetry.
% This function uses methods for describing the whole distribution of ISIs and thus
% is a crude way to specify the ISI distribution.
%
% Input:
%       <isi_lst>    - Vector or matrix containing interval durations between spikes.
%       <DIM>        - Dimension of trials, use 1 if in collums and 2 if in rows.
%
% Output:
%       <CV>         - coefficent of variation of ISIs
%       <mISI>       - mean ISI
%       <SK>         - coefficent of the asymmetry of the interval distribution ([3] Shinomoto et al. 2002)
%
% Note:
%       This method is not robust against changes of firing rate and should be used and interpreted with care.
%
% References:
%
% [1] Softy WR & Koch C (1993).
%     The High Irregular Firing of Cortical Cells is Inconsistent with Temporal Integration of Random EPSPs.
%     J NeuroSci 13(1), 334-350
%
% [2] Holt GR, Softky WR, Koch C, and Douglas RJ (1996)
%     Comparison of Discharge Variability In Vitro and In Vivo in Cat Visual Neurons.
%     J NeuroPhys, 75(5), 1806-1814
%
% [5] Shinomoto S, Sakai Y, Ohno H (2002).
%     Recording site dependence of the neuronal spiking statistics.
%     Biosystems, 67(1-3), 259-63
%
%
% wolf zinke, 03.09.2006

if (exist('DIM','var') == 0)
     DIM = 1;
elseif(isempty(DIM) == 1)
     DIM = 1;
end

%% determine first mean ISI and CV
isi_vec = isi_lst(:);
isi_vec(find(finite(isi_vec) == 0)) = [];   % Get rid of NaN's

if(length(isi_vec) < 1)   % There is nothing to do
    mISI = NaN;
    CV   = NaN;
    SK   = NaN;
    return;
end

mISI = mean(isi_vec);

% coefficent of variation as usually defined (well, as short as it in matlab gets)
CV = std(isi_vec) / mISI;

isi_mc = isi_vec - mISI;  % calculate the difference of ISIs and the mean ISI

% coefficent of the asymmetry of the interval distribution ([4] Shinomoto et al. 2002)
SK = mean(isi_mc.^3) / (mean(isi_mc.^2))^1.5;

