function [eta_s] = visc_birnbaum (f, eta_melt)

%model = 'birnbaum';
smax = 0.56;
gmax = 0.82;
Bgas = 1.98; 
Bsol = 2.74;

eta_s = eta_melt * ((1-f(3,:)).^(-Bgas)) .* ((1-f(1,:)/smax).^(-Bsol));

eta_s(f(1,:)>0.8*smax) = nan;
eta_s(f(3,:)>gmax) = nan;


end