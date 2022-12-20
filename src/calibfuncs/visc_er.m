function [eta] = visc_er (f, eta_melt, fsmax)

% Einstein-Roscoe effective liquid shear viscosity
if nargin==2, fsmax = 0.56; end

B1   = 4.0;  % theoretical value = 2.5
eta  = eta_melt .* max(1e-20,1-f(1,:)./fsmax).^-(fsmax.*B1);

eta(f(2,:)<(fsmax-0.2)) = nan;

end