function [kmvp,fout] = gcrit_parmi (f)
% 
% [phigcrit] = gcrit_parmi (f)
% 
% calculates the critical gas fraction from the Parmigiani et al. (2017)
% lattice-boltzmann models. 
% the models are valid for solid fractions between [0.39, 0.7].

% original from Parmigiani 2017
phigcrit = 2.75*f(1,:).^3 - 2.79*f(1,:).^2 + 0.6345*f(1,:) + 0.0997;

% alternative fit from degruyter 2019
% phigcrit = 0.7495*f(1,:).^3 - 0.4268*f(1,:).^2 - 0.1626*f(1,:) + 0.1478;

fout     = [f(1,:); 1-f(1,:)-phigcrit; phigcrit];

k  = 1e-4*(-0.0534*f(1,:).^3 + 0.1083*f(1,:).^2 - 0.0747*f(1,:) + 0.0176 );
kr = -2.1778*f(1,:).^4 + 5.1511*f(1,:).^3 - 4.5199*f(1,:).^2 + 1.7385*f(1,:) - 0.2461;
kmvp = k.*kr;

fout(:, f(1,:)<0.4 | f(1,:)>0.7) = nan;
kmvp(:, f(1,:)<0.4 | f(1,:)>0.7) = nan;


end