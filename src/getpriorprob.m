
function [pm] = getpriorprob (x, xd, distrib)
% 
% [pm] = getpriorprob (x, xd, distrib)
% 
% get the prior probability of a model given a distribution
% 
% INPUTS
% x         vector of model parameters [Nvar x 1]
% xd        matrix of properties of distribution [see below, size Nvar x 2]
% distrib   type of distribution, either 'normal' or 'uniform'
% 
% OUTPUT
% pm        joint probability of the model parameters in log space [scalar]
% 
% YQW, 20 April 2022


Nvar = length(x);
pvar = zeros(Nvar);

switch distrib
    case 'normal'   
        % where xd = [mean, standard deviation]
        pvar = normpdf(x(:), xd(:,1), xd(:,2));
        
    case 'uniform'  
        % where xd = [lower, upper]
        pvar = double(x(:)>=xd(:,1) & x(:)<=xd(:,2));
end

% get joint probability IN LOG SPACE
pm = sum( log(pvar) );

end

