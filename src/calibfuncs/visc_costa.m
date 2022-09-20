function [eta, zeta] =  visc_costa (f, eta_melt, eta_solid)
%
% [eta, zeta] =  visc_costa (f, eta_melt, eta_solid)
% 
% calculates the viscosity as a function of phase fraction using the 
% Costa et al. 2009 model. But this is still a model, so we need to fit it
% to the pure phase viscosities. 
% If only eta_melt is specified, I use a pre-fitted value for the plg-dac-mvp case. 
% If eta_solid is also specified, this function automatically fits xi,
% the parameter that controls the viscosity at solid fraction = 1.
% 
% 

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

feff = f(1:2,:)./sum(f(1:2,:),1);

% Costa+ effective mixture shear and compaction viscosities  (Costa et al., 2009)
B1        = 4.0;  % theoretical value = 2.5
phistar   = 0.62;
gamma     = 3.25;
delta     = 24;
xi        = 0.00148; % default, for plg-dac-mvp where solid visc = 1e16 Pa s.

% calculate xi!
if nargin==3
    xi = run_xifit(1, xi, phistar, gamma, delta, B1, eta_melt, eta_solid);
end
% for olv-bas: 4e-5 

eta  = costa_eta(feff(1,:), xi, phistar, gamma, delta, B1, eta_melt);
zeta = min(HUGE,max(TINY, eta./feff(2,:) ));

end

function [xifit] = run_xifit (f, xi0, phistar, gamma, delta, B1, eta_melt, eta_solid)

fprintf(1,'\n\n Running fsolve for costa model to fit \n');
fprintf(1,'  xi parameter to solid viscosity of %.0e Pa s.\n', eta_solid);

xifit   = fsolve(@(xif) log10(costa_eta(f, xif, phistar, gamma, delta, B1, eta_melt)) - log10(eta_solid), xi0);

fprintf(1, '\n\nBest fit xi = %.2e.\n\n', xifit);

end

function [eta] = costa_eta (f, xi, phistar, gamma, delta, B1, eta_melt)

fscl= f./phistar;
h   = (1-xi).*erf(sqrt(pi)./(2.*(1-xi)).*fscl.*(1+fscl.^gamma));
eta = eta_melt .* (1+fscl.^delta) ./ (1-h).^(B1*phistar);

end