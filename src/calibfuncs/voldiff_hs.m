function [kphi, Cseg] = voldiff_hs (fc, fd, etac, d0d)
% 
% [kphi, Cseg] = voldiff_hs (fc, fd, eta_melt, d0)
% 
% volume diffusivity and settling coefficient from hindered stokes settling
% 
% INPUTS
% fc    continuous phase (i.e. melt)
% fd    dispersed phase (i.e. crystals or bubbles)
% etac  viscosity of continuous phase
% d0d   granular scale of dispersed phase

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

% normalise phase fractions into two phases
fnorm = fc + fd;
fc    = fc./fnorm;
fd    = fd./fnorm;

% hindered-Stokes effective volume diffusivity
kphi = d0d.^2.*fc.^5./etac;  % drop prefactor 10*2/9

% Hindered-Stokes settling coefficient
Cseg = min(HUGE,max(TINY,  d0d.^2.*fd.*fc.^5./etac ));  % drop prefactor 2/9

% set values outside high liquid fraction regime to be nans
kphi(fc<0.60 | fc>0.999) = nan;
Cseg(fc<0.60 | fc>0.999) = nan;


end