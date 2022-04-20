%
% use this code to calibrate permission functions A, B, C for a plg-dac-mvp
% three-phase system
% NB: using log (base e) to define probabilities throughout the code
%
% YQW, 28 Jan 2022

clear all; clc; close all;
Addpaths
load('ocean.mat');
rng(5);

% addpath to mcmc
addpath('../../paramest/MCMC/');
addpath('../../paramest/catmip/');

%%
PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+4;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
% original
% A0 = [  0.60, 0.25, 0.30;
%         0.20, 0.20, 0.20;
%         0.20, 0.20, 0.20; ];  % permission slopes
% B0 = [  0.30, 0.15, 0.55;
%         0.48, 0.02, 0.50;
%         0.80, 0.08, 0.12; ];  % permission step locations
% C0 = [  0.20, 0.20, 0.20;
%         0.60, 0.60, 0.12;
%         0.20, 0.25, 0.50; ];  % permission step widths

% trial and error
A0 = [  0.60, 0.25, 0.30;
        0.20, 0.20, 0.20;
        0.20, 0.20, 0.20; ];  % permission slopes
B0 = [  0.30, 0.15, 0.55;
        0.48, 0.02, 0.50;
        0.48, 0.50, 0.02; ];  % permission step locations
C0 = [  0.20, 0.20, 0.20;
        0.60, 0.60, 0.12;
        0.07, 0.02, 0.20; ];  % permission step widths

vname = [strcat({'$A_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'});
         strcat({'$B_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'});
         strcat({'$C_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'})];

Ai =  1:9 ;
Bi = 10:18;
Ci = 19:27;

xccep = log10([A0(:); B0(:); C0(:)]);

%% define the models to calibrate against

dir('../src/calibfuncs/*.m')

Nf = 100;       % number of points to generate data on
Ns = 50;        % number of points to subsample for inversion

% set up phase fractions
f0   = SetUp3PhsMatrix(Nf);             % three-phase variation
fgl0 = SetUp3PhsMatrix(Nf, 1, 0.001);   % along gas-liquid axis
fgs0 = SetUp3PhsMatrix(Nf, 2, 0.001);   % along gas-solid axis
fsl0 = SetUp3PhsMatrix(Nf, 3, 0.001);   % along solid-liquid axis

% three phase viscosities
[eta_tb         ] = visc_truby(f0, eta0(2));
[eta_tb,f_eta_tb] = samplevalidvalues(log10(eta_tb),f0,Ns,1);

% viscosity along solid liquid axis
[ eta_co,   zeta_co] = visc_costa(fsl0, eta0(2));
  eta_co(fsl0(1,:)<0.3 | fsl0(1,:)>0.8) = nan;
 zeta_co(fsl0(1,:)<0.8) = nan;
[ eta_co,  f_eta_co] = samplevalidvalues(log10( eta_co),fsl0,Ns,1);
[zeta_co, f_zeta_co] = samplevalidvalues(log10(zeta_co),fsl0,Ns,1);

% hindered Stokes on dispersed solid phase
[vold_hs_s,   seg_hs_s] = voldiff_hs(fsl0(2,:), fsl0(1,:), eta0(2), d0(1));
[vold_hs_s,f_vold_hs_s] = samplevalidvalues(log10(vold_hs_s),fsl0,Ns,0);
[ seg_hs_s, f_seg_hs_s] = samplevalidvalues(log10( seg_hs_s),fsl0,Ns,0);

% hindered Stokes on dispersed bubble phase
[vold_hs_g,   seg_hs_g] = voldiff_hs(fgl0(2,:), fgl0(3,:), eta0(2), d0(3));
[vold_hs_g,f_vold_hs_g] = samplevalidvalues(log10(vold_hs_g),fgl0,Ns,0);
[ seg_hs_g, f_seg_hs_g] = samplevalidvalues(log10( seg_hs_g),fgl0,Ns,0);

% kozeny carman on melt percolation
[segKC_l          ] = perm_kc(fsl0(1,:), fsl0(2,:), eta0(2), d0(1));
[segKC_l,f_segKC_l] = samplevalidvalues(log10(segKC_l),fsl0,Ns,0);

% kozeny carman on mvp percolation
[segKC_g          ] = perm_kc(fgs0(1,:), fgs0(3,:), eta0(3), d0(1));
[segKC_g,f_segKC_g] = samplevalidvalues(log10(segKC_g),fgs0,Ns,0);

% critical gas fraction from parmigiani
[kmvp, f_gcrit] = gcrit_parmi(fsl0);
[kmvp, f_gcrit] = samplevalidvalues(log10(kmvp/eta0(3)),f_gcrit,Ns,0);


%% cast data into form for model outputs - all possible calib data
% possible options for dcat:
% etamix    mixture viscosity
% voldmix   mixture volume diffusion
% comp1     compaction coefficient of phase 1
% segr1     segregation coefficient of phase 1
% segr2     segregation coefficient of phase 2
% segr3     segregation coefficient of phase 3
% mvpcrt    critical mvp fraction for channelization onset

ftot = [f_gcrit, f_eta_tb, f_eta_co, f_zeta_co, f_segKC_l, f_seg_hs_s, f_seg_hs_g, f_segKC_g];
data = [  kmvp ,   eta_tb,   segKC_l,   eta_co,   zeta_co,   seg_hs_s,   seg_hs_g,   segKC_g];
dcat = repelem({'segr3','etamix','etamix','comp1','segr2','segr1','segr3','segr3'},1,Ns);
sigm = ones(size(data));

% plot where data is defined on tern axes
plotdataonternaxis(ftot, dcat);

% SaveFigure('Figures/possiblecalibmodels');

%% choose subset of calib models, due to issues with some models

ftot = [f_gcrit, f_eta_tb, f_segKC_l, f_eta_co, f_zeta_co, f_seg_hs_s];
data = [  kmvp ,   eta_tb,   segKC_l,   eta_co,   zeta_co,   seg_hs_s];
dcat = repelem({'segr3','etamix','segr2','etamix','comp1','segr1'},1,Ns);
sigm = 100*reshape([1,1,1,1,1,1].*ones(Ns,1),1,[])';

% plot where data is defined on tern axes
plotdataonternaxis(ftot, dcat);

%% define prior - choose one

distrib  = 'uniform';
xd       = zeros(3*NPHS^2,2);
xd(Ai,1) = -2;
xd(Bi,1) = -2;
xd(Ci,1) = -2;

% distrib     = 'normal';
% xd(:,1)     = xccep;
% xd(:,2)     = 0.25*ones(size(xccep));

%% define and check functions

psf = @(Nit ) sampleprior(xd, Nit, distrib);
dhf = @(m   ) calccoeffs(m, ftot, eta0, d0, dcat);
lmf = @(m   ) likelihoodfrommodel(m, ftot, eta0, d0, data, dcat, sigm);

dhf(xccep); lmf(xccep); psf(10);
plotfittodata(ftot, data, dhf(xccep), dcat, sigm);

%% run catmip

[xout, Pout, dhat, RunTime, allmodels] = catmip(psf,lmf,xd,'Niter',1000000);
% [xout] = normaliseB(xout, vname, 'log');

%% plot outputs

% visualize outputs
[xMAP] = plotdistribs(xout, Pout, vname, xd, distrib);
[A_MAP, B_MAP, C_MAP] = permvec2mat(xMAP, 'log')

%%

dhatMAP = dhf(xMAP);
plotfittodata(ftot, data, dhatMAP, dcat, sigm)
% SaveFigure('Figures/pdmcalib_mcmc_datafit');

fout = SetUp3PhsMatrix(200);
[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(fout, eta0, d0, A_MAP, B_MAP, C_MAP);

ax = Plot3PhasePerm(fout, Xf, PHS);
axes(ax(1)); hold on;
htb = ternplot(f_eta_tb(1,:), f_eta_tb(2,:), f_eta_tb(3,:), 'r+');
htb.ZData = 10.*ones(size(htb.XData));
% SaveFigure('Figures/pdmcalib_mcmc_Xf');

Plot3PhaseCoeff(fout, cat(3,Kv,Kf), 'scl', 'log', 'PHS', PHS, 'cfname', {'K_v/\phi','K_\phi/\phi'});
% SaveFigure('Figures/pdmcalib_mcmc_kphi');

Plot3PhaseCoeff(fout, cat(3,fout.^2./Cv,fout.^2./Cf), 'scl', 'log', 'PHS', PHS, 'cfname', {'\phi^2/C_v','\phi^2/C_\phi'});
% SaveFigure('Figures/pdmcalib_mcmc_segcompcoeff');





