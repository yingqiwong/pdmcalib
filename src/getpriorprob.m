
function [px] = getpriorprob (x, xd, distrib)
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
% px        joint probability of the model parameters in log space [scalar]
% 
% YQW, 20 April 2022


Nvar = length(x);
pvar = zeros(Nvar);

switch distrib
    case 'normal'   
        % where xd = [mean, standard deviation]
        pvar = - ( (x(:) - xd(:,1))./xd(:,2) ).^2;
                
    case 'uniform'  
        % where xd = [lower, upper]
        pvar = log( double(x(:)>=xd(:,1) & x(:)<=xd(:,2)) );
        
    case 'bounded_normal'
        pvar = - ( (x(:) - xd(:,1))./xd(:,2) ).^2;
        pvar(x(:)<xd(:,3) | x(:)>xd(:,4)) = -1e6;
end

% get joint probability 
px = sum( pvar );

end

