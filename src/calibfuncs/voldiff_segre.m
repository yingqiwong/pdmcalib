function [kphi] = voldiff_segre (fc, fd, etac, d0c)
% 
% [kphi] = voldiff_segre (fc, fd, eta_melt, d0_melt)
% 
% volume diffusivity and settling coefficient from hindered stokes settling
% 
% INPUTS
% fc    continuous phase (i.e. melt)
% fd    dispersed phase (i.e. crystals or bubbles
% etac  viscosity of continuous phase
% d0c   granular scale of continuous (?) phase

% normalise phase fractions into two phases
fnorm = fc + fd;
fc    = fc./fnorm;
fd    = fd./fnorm;

% Segre+ effective volume diffusivity  (Segre et al., 2001)
S    = ( (fd.^4 - 4.*fd.^3 + 4.*fd.^2 + 4.*fd + 1)./(fd-1).^4 ).^-1;
kphi = 1.*d0c.^2./etac.*(1-fd./0.71).^2.*S.^(1/2);  
% drop geometric prefactor 11.4*2/9

% set values outside high liquid fraction regime to be nans
kphi(fc<0.60 | fc>0.999) = nan;
end