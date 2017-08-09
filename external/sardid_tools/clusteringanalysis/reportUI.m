function ui = reportUI
% function that reports ui, that is whether Octave or matlab is in use
LIC = license('inuse');
ui = LIC.feature;
