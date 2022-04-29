%
% use this code to calibrate permission functions A, B, C for a plg-dac-mvp
% three-phase system
% NB: using log (base e) to define probabilities throughout the code
%
% YQW, 28 Jan 2022

clear all; clc; close all;

% add required paths, specify output directory
Addpaths; addpath('../../paramest/catmip/');

%% specify some inputs

outdir  = '../out/';        % output directory
runID   = 'test';           % name of this inversion

distrib = 'normal';        % 'uniform' or 'normal' prior distribution

% catmip options
Niter   = 1e4;              % number of iterations per catmip temperature
Nstep   = 400;              % number of steps in MCMC in catmip
pllopt  = 8;                % whether to run in parallel, number of workers

%% material properties

PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+4;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% initial set of parameters that I fit by trial and error
A0 = [  0.60, 0.25, 0.30;
        0.20, 0.20, 0.20;
        0.20, 0.20, 0.20; ];  % permission slopes
B0 = [  0.30, 0.15, 0.55;
        0.48, 0.02, 0.50;
        0.48, 0.50, 0.02; ];  % permission step locations
C0 = [  0.20, 0.20, 0.20;
        0.60, 0.60, 0.12;
        0.07, 0.02, 0.20; ];  % permission step widths

xccep = log10([A0(:); B0(:); C0(:)]);

%% specify variable indices

vname = [strcat({'$A_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'});
         strcat({'$B_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'});
         strcat({'$C_{'}, repmat(num2str((1:3)'),3,1),  repelem(num2str((1:3)'),3,1), {'}$'})];

Ai =  1:9 ;
Bi = 10:18;
Ci = 19:27;

%% define parameters of the prior distribution

switch distrib
    case 'uniform'
        xd(Ai,:) = [-2, 0].*ones(NPHS^2,1);
        xd(Bi,:) = [-2, 0].*ones(NPHS^2,1);
        xd(Ci,:) = [-2, 1].*ones(NPHS^2,1);
    case 'normal'
        xd(:,1)  = xccep;
        xd(:,2)  = 0.2*ones(size(xccep));
end

%% some other initialization for this inversion run

runID  = [runID '_' distrib];
rundir = [outdir runID '/'];

if ~exist(rundir, 'dir'), mkdir(rundir); end

% make a log file
logfile = [rundir runID '.log'];
if exist(logfile,'file'); delete(logfile); end
diary(logfile)

%% define the models to calibrate against

rng(5);
dir('../src/calibfuncs/*.m')

Nf = 100;       % number of points to generate data on
Ns = 40;        % number of points to subsample for inversion

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
  eta_co(fsl0(1,:)<0.3 | fsl0(1,:)>0.80) = nan;
 zeta_co(fsl0(1,:)<0.8 | fsl0(1,:)>0.99) = nan;
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


%% cast data into form for model outputs - choose some of the data from above
% possible options for dcat:
% etamix    mixture viscosity
% voldmix   mixture volume diffusion
% comp1     compaction coefficient of phase 1
% segr1     segregation coefficient of phase 1
% segr2     segregation coefficient of phase 2
% segr3     segregation coefficient of phase 3
% mvpcrt    critical mvp fraction for channelization onset

ftot = [f_gcrit, f_eta_tb, f_segKC_l, f_eta_co, f_zeta_co, f_seg_hs_s];
data = [  kmvp ,   eta_tb,   segKC_l,   eta_co,   zeta_co,   seg_hs_s];
dcat = repelem({'segr3','etamix','segr2','etamix','comp1','segr1'},1,Ns);
sigm = reshape([0.5,1,1,1,1,1].*ones(Ns,1),1,[])';

% plot where data is defined on tern axes
plotdataonternaxis(ftot, dcat);
SaveFigure([rundir runID '_data.pdf']);

%% define and check functions

psf = @(Nit) sampleprior(xd, Nit, distrib);
prf = @(xt ) getpriorprob(xt, xd, distrib);
dhf = @(xt ) calccoeffs(xt, ftot, eta0, d0, dcat);
lmf = @(xt ) likelihoodfrommodel(xt, ftot, eta0, d0, data, dcat, sigm);

dhf(xccep); lmf(xccep); psf(10);
plotfittodata(ftot, data, dhf(xccep), dcat, sigm);

%% run catmip

[xout, Pout, dhat, RunTime, allmodels] = catmip(...
    prf,psf,lmf,'Niter',Niter,'Nsteps',Nstep,'Parallel',logical(pllopt),'Ncores',pllopt);

%% plot outputs

% posterior distributions
[xMAP] = plotdistribs(xout, Pout, vname, xd, distrib);
SaveFigure([rundir runID '_posteriordistributions.pdf']);

% corner plots to see correlations between posterior distributions
figs = plotcorrelations(xout, Pout, vname, xd, distrib);
for fi = 1:length(figs)
    SaveFigure([rundir runID '_correlations_' num2str(fi) '.pdf'], figs(fi));
end

% print out MAP model
[A_MAP, B_MAP, C_MAP] = permvec2mat(xMAP, 'log')

% plot the fit to the datasets
plotfittodata(ftot, data, dhf(xMAP), dcat, sigm);
SaveFigure([rundir runID '_xMAP_fittodata.pdf']);

% now calculate connectivity and coefficients for plotting
fout = SetUp3PhsMatrix(200);
[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(fout, eta0, d0, A_MAP, B_MAP, C_MAP);

% plot permissions
ax = Plot3PhasePerm(fout, Xf, PHS);
SaveFigure([rundir runID '_xMAP_connectivity.pdf']);

% plot effective viscosity and volume diffusivity
Plot3PhaseCoeff(fout,cat(3,Kv,Kf),'scl','log','PHS',PHS,'cfname',{'(K_v/\phi)','(K_\phi/\phi)'});
SaveFigure([rundir runID '_xMAP_fluxcoeff.pdf']);

% plot segregation and compaction coefficients
Plot3PhaseCoeff(fout,cat(3,fout.^2./Cv,fout.^2./Cf),'scl','log','PHS',PHS,'cfname',{'(\phi^2/C_v)','(\phi^2/C_\phi)'});
SaveFigure([rundir runID '_xMAP_segcompcoeff.pdf']);


save([rundir runID '_par.mat'],'runID','distrib','Niter','Nstep','pllopt',...
    'rho0','eta0','d0','xccep','xd','ftot','data','dcat','sigm','xout','Pout');
diary off


