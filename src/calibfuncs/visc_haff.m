function [etafrac] = visc_haff (f)

% viscosity of a granular gas-like material from the Haff model as
% described in Andreotti et al. (2013) Chapter 5.2
% 
% note that the output is a fraction of viscosity normalised by the
% viscosity at a small particle fraction, because there are many (constant)
% parameters in the constitutive equations that we do not know.
% basically just compares the dependence on phase fractions.

fc = 0.64;
f0 = 0.1;

feff = f([1,3],:)./sum(f([1,3],:),1);

r0 = ( (fc./f0       ).^(1/3) - 1 ).^(-1);
r  = ( (fc./feff(1,:)).^(1/3) - 1 ).^(-1);

etafrac = 1 + r./r0;
etafrac(feff(1,:)>0.9*fc) = nan;

end