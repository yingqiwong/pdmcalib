function [eta, zeta] = visc_rudge (f, eta_s)

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

feff = f(1:2,:)./sum(f(1:2,:),1);

% Rudge effective solid shear and compaction viscosities  (Rudge, 2018)
nu       = min(HUGE,-1./log(feff(2,:)));
a1       = -4.28207265;
a2       = 7.36988663;
a3       = -4.98396638;
b1       = 0.97736898;
b2       = -1.76154195;
b3       = 2.63720462;


eta  = min(HUGE,max(TINY,  eta_s.*max(1e-20,1 + a1.*nu + a2.*nu.^2 + a3.*nu.^3)  ));
zeta = min(HUGE,max(TINY,  eta.*(160.*sqrt(2)./(139.*pi.*nu) + b1 + b2.*nu + b3.*nu.^2)  ));


% set phase fraction limits for this model
lmax = 0.30;
smax = 0.9999999;

eta( feff(2,:)>lmax | feff(1,:)>smax) = nan;
zeta(feff(2,:)>lmax | feff(1,:)>smax) = nan;


end
