function [kphi, Cseg] = voldiff_hs (f, eta_melt, d0)

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

% hindered-Stokes effective volume diffusivity
kphi = d0(2).^2.*f(2,:).^5./eta_melt;  % drop prefactor 10*2/9

% Hindered-Stokes settling coefficient
Cseg = min(HUGE,max(TINY,  d0(1).^2.*f(1,:).*f(2,:).^5./eta_melt ));  % drop prefactor 2/9


end