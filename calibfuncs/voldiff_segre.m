function [kphi] = voldiff_segre (f, eta_melt, d0_melt)


% Segre+ effective volume diffusivity  (Segre et al., 2001)
S    = ((f(1,:).^4 - 4.*f(1,:).^3 + 4.*f(1,:).^2 + 4.*f(1,:) + 1)./(f(1,:)-1).^4).^-1;
kphi = 1.*d0_melt.^2./eta_melt.*(1-f(1,:)./0.71).^2.*S.^(1/2);  
% drop geometric prefactor 11.4*2/9


end