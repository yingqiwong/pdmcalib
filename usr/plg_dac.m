%
% use this code to calibrate permission functions A, B, C for a two-phase
% system based on olv-bas mixture
% NB: using log (base e) to define probabilities throughout the code
%
% YQW, 7 Nov 2022

clear all; clc; close all;

% add required paths, specify output directory
Addpaths; addpath('../../paramest/catmip/');

%% specify some inputs

outdir  = '../out/';        % output directory
runID   = 'plgdac_largeNiter_2';    % name of this calibration

distrib = 'uniform';        % 'uniform' or 'normal' prior distribution

% catmip options
Niter   = 8^5;              % number of iterations per catmip temperature
Nstep   = 1e3;              % number of steps in MCMC in catmip
pllopt  = 8;                % whether to run in parallel, number of workers

%% material properties

PHS  = {'plg', 'dac'}; % phase names
NPHS = length(PHS);

rho0 = [ 2600; 2300];       % pure-phase densities
eta0 = [1e+16;1e+05];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents
costa_phistar = 0.52;

% set permission weight parameters for coefficient closure model
% define olv-bas model for comparison (not used in model)
A0 = [ 0.5989, 0.1772; 0.0397, 0.1182 ];  % permission slopes
B0 = [ 0.6870, 0.3130; 0.9998, 0.0002 ];  % permission step locations
C0 = [ 9.0105, 0.1592; 0.7249, 3.5524 ];  % permission step widths

xccep = log10([A0(:); B0(:); C0(:)]);

%% specify variable indices

vname = [strcat({'$A_{'}, repmat(num2str((1:NPHS)'),NPHS,1),  repelem(num2str((1:NPHS)'),NPHS,1), {'}$'});
         strcat({'$B_{'}, repmat(num2str((1:NPHS)'),NPHS,1),  repelem(num2str((1:NPHS)'),NPHS,1), {'}$'});
         strcat({'$C_{'}, repmat(num2str((1:NPHS)'),NPHS,1),  repelem(num2str((1:NPHS)'),NPHS,1), {'}$'})];

Ai = (1:NPHS^2)           ;
Bi = (1:NPHS^2) + 1*NPHS^2;
Ci = (1:NPHS^2) + 2*NPHS^2;

%% define parameters of the prior distribution

switch distrib
    case 'uniform'
        xd(Ai,:) = [-2, 0].*ones(NPHS^2,1);
        xd(Bi,:) = [-4, 0].*ones(NPHS^2,1);
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

rng(10);
dir('../src/calibfuncs/*.m')

Nf = 10001;     % number of points to generate data on
Ns = 200;        % number of points to subsample for inversion

fs  = linspace(0,1,Nf);
fsl = [fs; 1-fs];

% viscosity along solid liquid axis
[ eta_co,   zeta_co] = visc_costa(fsl, eta0(2), eta0(1), 'phistar', costa_phistar);
[ eta_co,  f_eta_co] = samplevalidvalues(log10( eta_co),fsl,Ns,0);
[zeta_co, f_zeta_co] = samplevalidvalues(log10(zeta_co),fsl,Ns,0);

% hindered Stokes on dispersed solid phase
[vold_hs,   seg_hs] = voldiff_hs(fsl(2,:), fsl(1,:), eta0(2), d0(1), costa_phistar-0.05);
[vold_hs,f_vold_hs] = samplevalidvalues(log10(vold_hs),fsl,Ns,0);
[ seg_hs, f_seg_hs] = samplevalidvalues(log10( seg_hs),fsl,Ns,0);

% Segre model for dispersed solid phase
[vold_sg          ] = voldiff_segre(fsl(2,:), fsl(1,:), eta0(2), d0(2), costa_phistar-0.05);
[vold_sg,f_vold_sg] = samplevalidvalues(log10(vold_sg),fsl,Ns,0);

% kozeny carman on melt percolation
[segKC        ] = perm_kc(fsl(1,:), fsl(2,:), eta0(2), d0(1));
[segKC,f_segKC] = samplevalidvalues(log10(segKC),fsl,Ns,0);


%% cast data into form for model outputs - choose some of the data from above
% possible options for dcat:
% etamix    mixture viscosity
% voldmix   mixture volume diffusion
% comp1     compaction  coefficient of phase 1
% segr1     segregation coefficient of phase 1
% segr2     segregation coefficient of phase 2
% segr3     segregation coefficient of phase 3
% mvpcrt    critical mvp fraction for channelization onset

ftot = [f_segKC, f_eta_co, f_seg_hs, f_vold_hs, f_vold_sg];
data = [  segKC,   eta_co,   seg_hs,   vold_hs,   vold_sg];
dcat = repelem({'segr2','etamix','segr1','voldmix','voldmix'},1,Ns);
sigm = 0.5*reshape([1,1,1,1,1].*ones(Ns,1),1,[])';

% plot where data is defined on tern axes
plotdataonlinaxis(ftot, dcat);
SaveFigure([rundir runID '_data.pdf']);

%% define and check functions

psf = @(Nit) sampleprior(xd, Nit, distrib);
prf = @(xt ) getpriorprob(xt, xd, distrib);
dhf = @(xt ) calccoeffs(xt, ftot, eta0, d0, dcat);
lmf = @(xt ) likelihoodfrommodel(xt, ftot, eta0, d0, data, dcat, sigm);

dhf(xccep); lmf(xccep); psf(10);
plotfittodata(ftot, Ns, data, dhf(xccep), dcat, sigm);

%% run catmip

[xout, Pout, dhat, RunTime, allmodels] = catmip(...
    prf,psf,lmf,'Niter',Niter,'Nsteps',Nstep,'Parallel',logical(pllopt),'Ncores',pllopt);

%% plot outputs

run('../src/output');

