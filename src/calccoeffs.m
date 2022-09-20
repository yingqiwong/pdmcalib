
function [dhat] = calccoeffs (x, f, eta0, d0, dcat)
% 
% [dhat] = calccoeffs (x, f, eta0, d0, dcat)
% 
% calculates the necessary coefficients in the forward model to compare
% with the calibration models
% 
% INPUTS
% x         model parameter vector [Nvar x 1]
% f         phase fractions at which to calculate coeffs 
% eta0      pure phase viscosity [NPHS x 1]
% d0        pure phase granular scale [NPHS x 1]
% dcat      category of models in data [Ndata x 1]
% 
% OUTPUT
% dhat      predicted data in the same form as calibration models [Ndata x 1]
% 
% YQW, 20 April 2022
% 

NPHS = size(f,1);
[A, B, C] = permvec2mat(x, 'log');

% set limit values
TINY = 1e-32;
HUGE = 1e+32;

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% get permission weights
F  = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
Sf = (F./B).^(1./C);  Sf = Sf./sum(Sf,2);
Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));
thtf = squeeze(prod(Mf.^Xf,2));

% get momentum and volume flux and transfer coefficients
Kv =    f .*eta0       .*thtv;
Kf =    f .*d0.^2./eta0.*thtf;
Cv = (1-f)./d0.^2.*Kv;
Cf = (1-f)./d0.^2.*Kf;

% control extreme values of coefficients
Kv = min(HUGE,max(TINY, Kv));
Kf = min(HUGE,max(TINY, Kf));
Cv = min(HUGE,max(TINY, Cv));
Cf = min(HUGE,max(TINY, Cf));

% calculate data from coeffs
eta_mix = sum(Kv,1);    % mixture viscosity
voldiff = sum(Kf,1);    % mixture volume diffusion coefficient
segcoef = f.^2./Cv;     % phase-wise segregation coefficient
cmpcoef = f.^2./Cf;     % phase-wise compaction coefficient

% assemble data vector
dhat = nan(1,size(f,2));
dhat(strcmp(dcat,'etamix' )) = log10( eta_mix(1,strcmp(dcat,'etamix' )) );
dhat(strcmp(dcat,'voldmix')) = log10( voldiff(1,strcmp(dcat,'voldmix')) );
dhat(strcmp(dcat,'comp1'  )) = log10( cmpcoef(1,strcmp(dcat,'comp1'  )) );
dhat(strcmp(dcat,'segr1'  )) = log10( segcoef(1,strcmp(dcat,'segr1'  )) );
dhat(strcmp(dcat,'segr2'  )) = log10( segcoef(2,strcmp(dcat,'segr2'  )) );

if NPHS>2
dhat(strcmp(dcat,'segr3'  )) = log10( segcoef(3,strcmp(dcat,'segr3'  )) );
end
end
