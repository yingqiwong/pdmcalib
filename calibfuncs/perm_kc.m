function [k] = perm_kc (f, eta, d0)

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

% Kozeny-Carman Darcy percolation coefficient
k = min(HUGE,max(TINY,  d0.^2./50.*(f(2,:)-0.001).^2.75.*f(1,:).^-2./eta  ));



end