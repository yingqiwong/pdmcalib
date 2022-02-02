function [eta] = visc_er (f, eta_melt)

% Einstein-Roscoe effective liquid shear viscosity
h    = 0.56;
B1   = 4.0;  % theoretical value = 2.5
eta  = eta_melt .* max(1e-20,1-f(1,:)./h).^-(h.*B1);

end