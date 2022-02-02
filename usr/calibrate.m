%
% use this code to calibrate permission functions A, B, C for a plg-dac-mvp
% three-phase system
% NB: using log (base e) to define probabilities throughout the code
% 
% YQW, 28 Jan 2022

clear all; clc;
Addpaths

% addpath to mcmc
addpath('../../../../MATLAB/lib/paramest/MCMC/');

%%
PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
% original
A0 = [  0.60, 0.25, 0.30;
        0.20, 0.20, 0.20;
        0.20, 0.20, 0.20; ];  % permission slopes
B0 = [  0.30, 0.15, 0.55;
        0.48, 0.02, 0.50;
        0.80, 0.08, 0.12; ];  % permission step locations
C0 = [  0.20, 0.20, 0.20;
        0.60, 0.60, 0.12;
        0.20, 0.25, 0.50; ];  % permission step widths

vname = [strcat({'$A_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'});
         strcat({'$B_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'});
         strcat({'$C_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'})];

%% define prior

Nm    = 3*NPHS^2;
mbnds = repmat([0,1], Nm, 1);
mstep = 0.05*diff(mbnds,[],2); 

m0    = mbnds(:,1) + diff(mbnds,[],2).*rand(Nm,1);
% m0    = [A0(:); B0(:); C0(:)];

[A, B, C] = permvec2mat(m0, NPHS);

%% define the data

Npm = 20;

f0 = SetUp3PhsMatrix(200);
f0(:,isnan(f0(1,:))) = [];
f0(:,f0(1,:)==1) = [];      f0(:,f0(1,:)==0) = [];
f0(:,f0(2,:)==1) = [];      f0(:,f0(2,:)==0) = [];
f0(:,f0(3,:)==1) = [];      f0(:,f0(3,:)==0) = [];


% viscosity at liquid-rich end
eta_bb = visc_birnbaum(f0, eta0(2));
[eta_bb,fb] = samplevalidvalues(log10(eta_bb),f0,Npm);

% viscosity at solid-rich end
% can't do this because you will have multi-valued solutions for different
% mvp. need to specify that mvp = one constant small value and fit the
% eta_r, zeta_r to that
[eta_r, zeta_r]  = visc_rudge(f0, eta0(1));
[eta_r ,fr_eta ] = samplevalidvalues(log10( eta_r),f0,Npm);
[zeta_r,fr_zeta] = samplevalidvalues(log10(zeta_r),f0,Npm);

% critical gas fraction from parmigiani
phigcrit = gcrit_parmi(f0);
[phigcrit,fp] = samplevalidvalues(phigcrit,f0,Npm);
fp(3,:) = phigcrit; fp(2,:) = 1 - fp(1,:) - fp(3,:);

ftot = [fb    ,fr_eta,fr_zeta,              fp];
data = [eta_bb, eta_r, zeta_r, 0.1*ones(1,Npm)];
dcat = [ones(1,Npm), 2*ones(1,Npm), 3*ones(1,Npm), 4*ones(1,Npm)];
sigma = ones(size(data));


%% define and check functions

prf = @(m) prior(m,mbnds);
lkf = @(dhat) likelihood(dhat, data, sigma);
dhf = @(m) calcdhat(m, ftot, eta0, d0, dcat);

prf(m0);
dhattmp = dhf(m0);
lkf(dhattmp);


%%

Niter = 1e5;
[x_keep,P_keep,count] = mcmc(dhf,prf,lkf,m0,mstep,mbnds,Niter);
PlotMCMCAnalytics(x_keep,P_keep,mbnds,count,0.2*Niter,vname);

%% functions

function [A, B, C] = permvec2mat (m, NPHS)
% manipulate model to permission curve functions
A = reshape(m(          1:NPHS^2  ), NPHS, NPHS);
B = reshape(m(1*NPHS^2+(1:NPHS^2) ), NPHS, NPHS); B = B./sum(B,2);
C = reshape(m(2*NPHS^2+(1:NPHS^2) ), NPHS, NPHS);
end

function [f, data] = setupdata (eta0, d0, Npm)

f0 = SetUp3PhsMatrix(100);
f0(:,isnan(f0(1,:))) = [];

% viscosity at liquid-rich end
[eta_bb,fb] = visc_birnbaum(f0, eta0(2));
fb = fb(:,randi(length(eta_bb),Npm,1));
eta_bb = eta_bb(randi(length(eta_bb),Npm,1));

% viscosity at solid-rich end
[eta_r, zeta_r] = visc_rudge(f0, eta0(1));
fr = f0;
fr(  :,isnan(eta_r)) = []; 
eta_bb(isnan(eta_bb)) = [];
fb = fb(:,randi(length(eta_bb),Npm,1));
eta_bb = eta_bb(randi(length(eta_bb),Npm,1));

end

function [vout, fout] = samplevalidvalues (v, f, Npm)

f(:,isnan(v(1,:))) = [];
v(:,isnan(v(1,:))) = [];

ind  = randperm(size(v,2), Npm);
fout = f(:,ind);
vout = v(:,ind);

end

function [pm] = prior (m, mbnds)
% define prior probability
pm = sum(log( double(m>=mbnds(:,1) & m<=mbnds(:,2) )) );
end

function ldm = likelihood (dhat, data, sigma)
% calc data likelihood (how well the predicted A,B,C match constraints)
ldm = sum(- 0.5*((data(:) - dhat(:))./sigma(:)).^2 );
end

function [dhat] = calcdhat (m, f, eta0, d0, dcat)
% calc coeffs and shape into dhat vector

NPHS = size(f,1);
[A, B, C] = permvec2mat(m, NPHS);

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

dhat          = nan(1,size(f,2));
dhat(dcat==1) = log10( sum(Kv(:,dcat==1), 1) );
dhat(dcat==2) = log10( sum(Kv(:,dcat==2), 1) );
dhat(dcat==3) = log10( sum( f(:,dcat==3).^2./Cf(:,dcat==3), 1) );
dhat(dcat==4) = squeeze(Xf(3,3,dcat==4));

% dhat(isnan(dhat)) = 1e10;

end



