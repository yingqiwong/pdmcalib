function [eta_s] = visc_truby (f, eta_melt)


% model = 'truby';
smax = 0.593;
gmax = 0.3;
Bgas = 1; 
Bsol = 2;

eta_s = eta_melt * ((1-f(3,:)).^(-Bgas)) .* ((1-f(1,:)/smax).^(-Bsol));

eta_s(f(1,:)>0.6*smax) = nan;
eta_s(f(3,:)>0.8*gmax) = nan;


end