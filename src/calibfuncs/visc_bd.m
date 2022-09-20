function [eta] = visc_bd (f, eta_melt)
%
% [eta] = visc_bd (f, eta_melt)
% 
% calculates the relative viscosity of a melt containing bubbles from
% experiments on porous glasses of Bagdassarov and Dingwell (1992)
%    


feff = f;
feff(2:3,:) = f(2:3,:)./sum(f(2:3,:),1);

b     = 22.4;
eta_b = 1 ./ (1 + b*feff(3,:));
eta   = eta_b.*eta_melt;


end