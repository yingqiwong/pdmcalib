function [eta] = visc_mader (f, eta_melt, Ca)
%
% [eta] = visc_mader (f, eta_melt, Ca)
% 
% calculates the relative viscosity of a melt containing bubbles from
% from theory of Llewellin et al. (2002) and extended by Mader et al.
% (2013)
%    


feff = f;
feff(2:3,:) = f(2:3,:)./sum(f(2:3,:),1);

K = 6/5;
m = 2;

eta_b_0   = (1 - feff(3,:)).^(-1);
eta_b_inf = (1 - feff(3,:)).^( 5/3);

eta_b = eta_b_inf + (eta_b_0 - eta_b_inf)./( 1 + (K*Ca)^m );
eta   = eta_b.*eta_melt;


end