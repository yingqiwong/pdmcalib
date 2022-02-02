function [phigcrit] = gcrit_parmi (f)
% 
% [phigcrit] = gcrit_parmi (f)
% 
% calculates the critical gas fraction from the Parmigiani et al. (2017)
% lattice-boltzmann models. 
% the models are valid for solid fractions between [0.39, 0.7].

phigcrit = 0.7495*f(1,:).^3 - 0.4268*f(1,:).^2 - 0.1626*f(1,:) + 0.1478;
phigcrit(f(1,:)<0.39 | f(1,:)>0.7) = nan;


end