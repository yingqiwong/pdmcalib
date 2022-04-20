
function [A, B, C] = permvec2mat (x, scale)
% 
% [A, B, C] = permvec2mat (x)
% 
% manipulate vector of model fitting parameters to the usual permission
% matrices A, B, C
% 
% INPUTS
% x         model vector corresponding to [A(:); B(:); C(:)]
%           size [(3 NPHS^2) x 1]
% 
% OUTPUTS
% A, B, C   fitting parameters for permission weights [NPHS x NPHS each]
% 
% YQW, 20 April 2022

if nargin==1, scale = []; end   % default is linear scale

NPHS = sqrt(length(x)/3);

A = reshape(x(          1:NPHS^2 ), NPHS, NPHS);    % slope
B = reshape(x(1*NPHS^2+(1:NPHS^2)), NPHS, NPHS);    % step location
C = reshape(x(2*NPHS^2+(1:NPHS^2)), NPHS, NPHS);    % step width

if strcmp(scale, 'log')
    A = 10.^A; 
    B = 10.^B;
    C = 10.^C;
end

% normalize B (sum across each row = 1)
B = B./sum(B,2);

end
