
function [xout] = normaliseB (xin, vname, scale)
%
% [xout] = normaliseB (xin, vname)
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

NPHS = sqrt(length(vname)/3);

% collect B's
Bvec = xin(:,contains(vname, 'B'));
if strcmp(scale, 'log'), Bvec = 10.^B; end

% now normalize B across each row to be 1. Loop over rows
for pi = 1:NPHS
    ri = pi + [0,3,6];
    Bvec(:,ri) = Bvec(:,ri)./sum(Bvec(:,ri),2);
end

% prepare output
xout = xin;
if strcmp(scale, 'log')
    xout(:,contains(vname, 'B')) = log10(Bvec);
else
    xout(:,contains(vname, 'B')) = Bvec;
end
end