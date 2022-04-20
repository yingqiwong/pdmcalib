
function ldm = likelihood (dhat, data, sigma)
% 
% ldm = likelihood (dhat, data, sigma)
% 
% calc data likelihood (how well the predicted A,B,C match constraints)
% given the predicted data from the model and the models against which to
% calibrate to
% 
% INPUTS
% dhat      predicted data from permission model 
% data      parameterization of models for calibration
% sigma     how well we think we should fit the models 
%   ^^^ all three input matrices can be of any size but MUST MATCH
% 
% OUTPUT
% ldm       likelihood of the set of model parameters IN LOG SPACE,
%           assuming gaussian errors
% 
% YQW, 20 April 2022
% 

ldm = sum(-((data(:) - dhat(:))./sigma(:)).^2 );

end
