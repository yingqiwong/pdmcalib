function [eta, zeta] =  visc_costa (f, eta_melt, eta_solid, varargin)
%
% [eta, zeta] =  visc_costa (f, eta_melt, eta_solid, disagg)
% 
% calculates the viscosity as a function of phase fraction using the 
% Costa et al. 2009 model. But this is still a model, so we need to fit it
% to the pure phase viscosities. 
% If only eta_melt is specified, I use a pre-fitted value for the plg-dac-mvp case. 
% If eta_solid is also specified, this function automatically fits xi,
% the parameter that controls the viscosity at solid fraction = 1.
% 

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

% normalize solid (1) and melt (2) fractions
feff = f(1:2,:)./sum(f(1:2,:),1);

% load parameters
[B1, phistar, gamma, delta, xi, fitxi] = defaultparams(varargin{:});

% calculate xi! for olv-bas: 4e-5
if (fitxi)
    fprintf(1,'\n\n Running fsolve for costa model to fit \n');
    fprintf(1,'  xi to solid viscosity of %.0e Pa s.\n', eta_solid);
    xi = fsolve(@(xif) log10(costa_eta(1, eta_melt, B1, phistar, gamma, delta, xif)) - log10(eta_solid), xi);
    fprintf(1, '\n\nBest fit xi = %.2e.\n\n', xi);
end

eta  = costa_eta(feff(1,:), eta_melt, B1, phistar, gamma, delta, xi);
zeta = min(HUGE,max(TINY, eta./feff(2,:) ));

end

function [B1, phistar, gamma, delta, xi, fitxi] = defaultparams (varargin)

p = inputParser;
p.addParameter('B1'     ,  4.0     , @isnumeric);
p.addParameter('phistar',  0.62    , @isnumeric);
p.addParameter('gamma'  ,  3.25    , @isnumeric);
p.addParameter('delta'  ,  24.0    , @isnumeric);
p.addParameter('xi'     ,  4e-5    , @isnumeric);
p.addParameter('fitxi'  ,  true    , @islogical);
p.addParameter('fitdelta', true    , @islogical);
p.parse(varargin{:}); 

B1      = p.Results.B1;
phistar = p.Results.phistar;
gamma   = p.Results.gamma;
delta   = p.Results.delta;
xi      = p.Results.xi;
fitxi   = p.Results.fitxi;

end

function [eta] = costa_eta (f, eta_melt, B1, phistar, gamma, delta, xi)

fscl= f./phistar;
h   = (1-xi).*erf(sqrt(pi)./(2.*(1-xi)).*fscl.*(1+fscl.^gamma));
eta = eta_melt .* (1+fscl.^delta) ./ (1-h).^(B1*phistar);

end
