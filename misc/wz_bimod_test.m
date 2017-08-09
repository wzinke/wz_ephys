function [b] = wz_bimod_test(x)
%
%  code from:
%  JB Freeman & R Dale
%  Assessing bimodality to detect the presence of a dual cognitive process
%  Behav Res Meth 2012
%  doi: 10.3758/s13428-012-0225-x
%
%  and:
%  R Pfister, KA Schwarz, M Janczyk, R Dale and JB Freeman
%  Good things peak in pairs: a note on the bimodality coefficient
%  Front. Psychol., 02 October 2013
%  doi: 10.3389/fpsyg.2013.00700
%

%m3 = skew
%m4 = kurt
%n = data size
m3 = skewness(x);
m4 = kurtosis(x, 0) - 3;
n = length(x);
b=(m3^2+1) / (m4 + 3 * ( (n-1)^2 / ((n-2)*(n-3)) ));