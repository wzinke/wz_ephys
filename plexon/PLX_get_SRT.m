function EY = PLX_get_SRT(EyeX, EyeY)
% Get a better estimate of a saccadic response.
%
% DESCRIPTION
%   The onset of a saccadic response is identified by smoothing the eyw traces
%   with a loess filter, identifying the time when a sufficient fast and large
%   eye movement exceeds a threshold (common approach) and then lok reversely
%   when this movement starts.
%
%
% SYNTAX
%
%
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 20-Feb-2015 by wolf zinke
