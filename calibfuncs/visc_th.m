function [eta, zeta, XF] = visc_th (f, eta_s)

% Takei-Holtzman connectivity  (Takei & Holtzmann, 2009)
a  = 2.3;
b  = 0.5;
xf = max(1e-32,1-a.*f(2,:).^b);
XF = log(1e18*xf.^2)./log(1e18);

% Takey-Holtzman effective solid shear and compaction viscosities  (Takei & Holtzmann, 2009)
c    = 0.2;
d    = 2;
eta  = c.*eta_s.*xf.^d;
zeta = c.*eta_s.*xf.^d;



end