function [eta] = visc_dr (f, eta_melt)
%
% [eta] = visc_dr (f, eta_melt)
% 
% calculates the relative viscosity of a melt containing bubbles from
% experiments on porous glasses of Ducamp and Raj (1989)
%    

feff = f;
feff(2:3,:) = f(2:3,:)./sum(f(2:3,:),1);


b     = 3;
eta_b = exp( -b*feff(3,:)./(1-feff(3,:)) );
eta   = eta_b.*eta_melt;


end