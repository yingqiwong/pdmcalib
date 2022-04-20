function [eta, zeta] =  visc_costa (f, eta_melt)

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

feff = f(1:2,:)./sum(f(1:2,:),1);

% Costa+ effective mixture shear and compaction viscosities  (Costa et al., 2009)
B1        = 4.0;  % theoretical value = 2.5
phistar   = 0.62;
phiscaled = feff(1,:)./phistar;
gamma     = 3.25;
delta     = 24;
xi        = 0.00148; %4e-5 in olv-bas;  % or 8e-8 (for B1 = 2.5)
h         = (1-xi).*erf(sqrt(pi)./(2.*(1-xi)).*phiscaled.*(1+phiscaled.^gamma));


eta  = eta_melt .* (1+phiscaled.^delta) ./ (1-h).^(B1*phistar);
zeta = min(HUGE,max(TINY, eta./feff(2,:) ));

end