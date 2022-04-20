function [eta] = visc_hk (f, eta_s)


% Hirth-Kohlstedt effective solid shear viscosity  (Hirth & Kohlstedt, 2003)
lam  = 27;
eta  = eta_s.*exp(-lam.*f(2,:));

end