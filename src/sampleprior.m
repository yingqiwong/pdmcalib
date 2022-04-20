
function [ms] = sampleprior (xd, Niter, distrib)
% 
% [ms] = sampleprior (xd, Niter, distrib)
% 
% obtain samples from a uniformly or normally distributed prior
% 
% INPUTS
% xd        matrix of properties of distribution [see below, size Nvar x 2]
% Niter     number of samples to return [scalar]
% distrib   type of distribution, either 'normal' or 'uniform'
% 
% OUTPUT
% ms        matrix of model parameters [Niter x Nvar]
% 
% YQW, 20 April 2022


% number of model parameters
Nvar = size(xd,1);

switch distrib
    
    case 'normal'
        % where xd = [mean, standard deviation]
        ms = xd(:,1) +  xd(:,2).*randn(Nvar,Niter);
        
    case 'uniform'
        % where xd = [lower, upper]
        ms = xd(:,1) + diff(xd,[],2).*rand(Nvar,Niter);
        
end

ms = ms';   % the output size has to be Niter x Nvar

end