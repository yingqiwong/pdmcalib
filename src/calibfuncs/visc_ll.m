function [eta] = visc_ll (f, eta_melt)
%
% [eta] = visc_ll (f, eta_melt)
% 
% calculates the relative viscosity of a melt containing bubbles from
% theoretical treatment of Llewellin et al. 2002
%    


feff = f;
feff(2:3,:) = f(2:3,:)./sum(f(2:3,:),1);

b     = 5/3;
eta_b = 1 - b*feff(3,:);
eta   = eta_b.*eta_melt;

eta(feff(3,:)>0.5) = nan;


end