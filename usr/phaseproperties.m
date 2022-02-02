


PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
% original
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% mess about with original
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.10, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.60, 0.20, 0.20; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% inspired by olv-bas
% A = [ 0.69, 0.18, 0.30; 0.54, 0.18, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
% B = [ 0.55, 0.30, 0.15; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
% C = [ 0.10, 0.18, 0.20; 0.82, 0.40, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% inspired by parmigiani 2017
% A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
% B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
% C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 1.20, 0.25, 0.50; ];  % permission step widths
