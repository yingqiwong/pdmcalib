function [eta] = visc_rahaman (f, eta_melt)
%
% [eta] = visc_rahaman (f, eta_melt)
% 
% calculates the relative viscosity of a melt containing bubbles from
% experiments on porous glasses of Rahaman et al. (1987)
%    

feff = f(2:3,:)./sum(f(2:3,:),1);


b     = 11.2;
eta_b = exp( -b*feff(3,:) );
eta   = eta_b.*eta_melt;


end