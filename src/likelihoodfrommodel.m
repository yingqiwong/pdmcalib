
function [ldm,dhat] = likelihoodfrommodel (x, f, eta0, d0, data, dcat, sigma)
% 
% [ldm,dhat] = likelihoodfrommodel (x, f, eta0, d0, data, dcat, sigma)
% 
% calc data likelihood (how well the predicted A,B,C match constraints)
% given the model parameters 
% basically combines the calccoeffs and likelihood functions, needed in
% this form to run catmip
% 
% INPUTS
% x         model parameter vector [Nvar x 1]
% f         phase fractions at which to calculate coeffs 
% eta0      pure phase viscosity [NPHS x 1]
% d0        pure phase granular scale [NPHS x 1]
% data      parameterization of models for calibration
% dcat      category of models in data
% sigma     how well we think we should fit the models 
%   ^^^ data, dcat and sigma can be of any size but MUST MATCH dhat
% 
% OUTPUT
% ldm       likelihood of the set of model parameters IN LOG SPACE,
%           assuming gaussian errors
% dhat      predicted data from calccoeffs
% 
% YQW, 20 April 2022
% 


dhat = calccoeffs(x, f, eta0, d0, dcat);     % run forward model
ldm  = likelihood(dhat, data, sigma);       % calculate likelihood

end
