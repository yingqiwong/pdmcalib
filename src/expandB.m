
function [xout] = expandB (xin, Ai, Bi, Ci, scale)
%
% [xout] = expandB (xin, vname)
%
% recalculate the permission step locations of a matrix of models so that
% they are normalized to 1, since the sum across each row of B should be 1
% but this is not possible to be enforced in the parameter estimation
% 
% INPUTS
% xin       matrix of models [Niter x Nvar]
%           where Nvar = 3 x NPHS^2
% vname     vector of parameter names [Nvar x 1]
% 
% OUTPUTS
% xout      output matrix of models where B values are normalized to 1 [Niter x Nvar]
% 
% YQW, 20 April 2022

if nargin < 3, scale = []; end

Nit  = size(xin,1);
NPHS = sqrt( length(Ai) );

% collect B's
Bvec = [xin(:,Bi), zeros(Nit,NPHS)];
if strcmp(scale, 'log'), Bvec = 10.^Bvec; end

% now normalize B across each row to be 1. Loop over rows
for iphs = 1:NPHS
    Bvec(:,2*NPHS+iphs) = 1 - Bvec(:,iphs) - Bvec(:,NPHS+iphs);
end


if strcmp(scale, 'log'), Bvec = log10(Bvec); end

% prepare output
xout = [xin(:,Ai), Bvec, xin(:,Ci)];

end