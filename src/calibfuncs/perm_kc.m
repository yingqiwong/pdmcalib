function [Cseg] = perm_kc (fmat, fperc, eta_perc, d0_mat)

% 
% [Cseg] = perm_kc (fmat, fperc, eta_perc, d0_mat)
% 
% segregation coefficient from kozeny-carman relationship
% 
% INPUTS
% fmat      matrix phase (i.e. solid)
% fperc     percolating phase (i.e. liquid or gas)
% eta_perc  viscosity of percolating phase
% d0_mat    granular scale of matrix phase


% set limit values
TINY = 1e-32;
HUGE = 1e+32;

% normalise phase fractions into two phases
ftot  = fmat + fperc;
fmat  = fmat./ftot;
fperc = fperc./ftot;

% Kozeny-Carman Darcy percolation coefficient
Cseg = min(HUGE,max(TINY,  d0_mat.^2./50.*(fperc-0.001).^2.75.*fmat.^-2./eta_perc  ));

% set values outside high solid fraction regime to be nans
Cseg(fmat<0.80 | fmat>0.97) = nan;



end