
function [x] = mat2permvec (A, B, C, scale)
% 
% [x] = mat2permvec (A, B, C, scale)
% 
% cast the permission weights fitting parameters into a vector
% 
% INPUTS
% A, B, C   fitting parameters for permission weights [NPHS x NPHS each]
% 
% OUTPUTS
% x         model vector corresponding to [A(:); B(:,1:(NPHS-1)); C(:)]
%           size [(NPHS^2 + (NPHS-1)xNPHS + NPHS^2) x 1]
% 
% YQW, 29 April 2022

if nargin==1, scale = []; end   % default is linear scale

NPHS = size(A,1);

x = [ A(:); 
      reshape(B(:,1:(NPHS-1)),[],1); 
      C(:)  ];

if strcmp(scale, 'log'), x = log10(x); end

end
