%
% use this code to calibrate permission functions A, B, C for a two-phase
% system based on olv-bas mixture
% NB: using log (base e) to define probabilities throughout the code
%
% YQW, 29 April 2022

clear all; clc; close all;

% add required paths, specify output directory
Addpaths; addpath('../../paramest/catmip/');

%% specify some inputs

outdir  = '../out/';        % output directory
runID   = 'olvbas_longchain';           % name of this calibration

distrib = 'uniform';        % 'uniform' or 'normal' prior distribution

% catmip options
Niter   = 1e4;              % number of iterations per catmip temperature
Nstep   = 1e3;              % number of steps in MCMC in catmip
pllopt  = 8;                % whether to run in parallel, number of workers

%% material properties

PHS  = {'olv', 'bas'}; % phase names
NPHS = length(PHS);

rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A0 = [ 0.6945, 0.1832; 0.5360, 0.1834;];  % permission slopes
B0 = [ 0.6906, 0.3094; 0.9993, 0.0007;];  % permission step locations
C0 = [ 0.6889, 0.1750; 0.8154, 1.5642;];  % permission step widths

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

rng(5);
dir('../src/calibfuncs/*.m')

Nf = 10001;     % number of points to generate data on
Ns = 400;        % number of points to subsample for inversion

fs  = linspace(0,1,Nf);
fsl = [fs; 1-fs];

% viscosity along solid liquid axis
[ eta_co,   zeta_co] = visc_costa(fsl, eta0(2), eta0(1));
[ eta_co,  f_eta_co] = samplevalidvalues(log10( eta_co),fsl,Ns,1);
[zeta_co, f_zeta_co] = samplevalidvalues(log10(zeta_co),fsl,Ns,1);

% hindered Stokes on dispersed solid phase
[vold_hs,   seg_hs] = voldiff_hs(fsl(2,:), fsl(1,:), eta0(2), d0(1));
[vold_hs,f_vold_hs] = samplevalidvalues(log10(vold_hs),fsl,Ns,0);
[ seg_hs, f_seg_hs] = samplevalidvalues(log10( seg_hs),fsl,Ns,0);

% Segre model for dispersed solid phase
[vold_sg          ] = voldiff_segre(fsl(2,:), fsl(1,:), eta0(2), d0(2));
[vold_sg,f_vold_sg] = samplevalidvalues(log10(vold_sg),fsl,Ns,0);

% kozeny carman on melt percolation
[segKC        ] = perm_kc(fsl(1,:), fsl(2,:), eta0(2), d0(1));
[segKC,f_segKC] = samplevalidvalues(log10(segKC),fsl,Ns,0);


%% cast data into form for model outputs - choose some of the data from above
% possible options for dcat:
% etamix    mixture viscosity
% voldmix   mixture volume diffusion
% comp1     compaction coefficient of phase 1
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

xnb = normaliseB(xout, vname, 'log');

% posterior distributions
[xMAP] = plotdistribs(xnb, Pout, vname, xd, distrib);
SaveFigure([rundir runID '_posteriordistributions.pdf']);

% corner plots to see correlations between posterior distributions
figs = plotcorrelations(xnb, Pout, vname, xd, distrib);
for fi = 1:length(figs)
    SaveFigure([rundir runID '_correlations_' num2str(fi) '.pdf'], figs(fi));
end

% print out MAP model
[A_MAP, B_MAP, C_MAP] = permvec2mat(xMAP, 'log')

% plot the fit to the datasets
plotfittodata(ftot, Ns, data, dhf(xccep), dcat, sigm);
SaveFigure([rundir runID '_xMAP_fittodata.pdf']);

%%  now calculate connectivity and coefficients for plotting

[dsc, KvMAP, KfMAP, CvMAP, CfMAP, XfMAP] = SegCompLength(fsl, eta0, d0, A_MAP, B_MAP, C_MAP);
[dsc, Kv0  , Kf0  , Cv0  , Cf0  , Xf0  ] = SegCompLength(fsl, eta0, d0, A0   , B0   , C0   );


figure;
hAx = setupaxes(3,2,'gaph',0.8,'gapw',0.3);
set(gcf,'defaultaxescolororder',repmat(lines(2),2,1));

axes(hAx(1));
plot(fsl(2,:), squeeze(XfMAP(1,1,:)), '-'); hold on;
plot(fsl(2,:), squeeze(XfMAP(2,2,:)), '-'); 
plot(fsl(2,:), squeeze(  Xf0(1,1,:)), ':');
plot(fsl(2,:), squeeze(  Xf0(2,2,:)), ':'); hold off;
ylabel('intra-phase weights');

axes(hAx(2));
plot(fsl(2,:), squeeze(XfMAP(1,2,:)), '-'); hold on;
plot(fsl(2,:), squeeze(XfMAP(2,1,:)), '-'); 
plot(fsl(2,:), squeeze(  Xf0(1,2,:)), ':'); 
plot(fsl(2,:), squeeze(  Xf0(2,1,:)), ':'); hold off;
ylabel('inter-phase weights');
set(gca,'YAxisLocation','right');

axes(hAx(3));
semilogy(fsl(2,:), KvMAP./fsl, '-'); hold on;
semilogy(fsl(2,:),   Kv0./fsl, ':'); hold off;
ylabel('effective viscosity [Pa s]');

axes(hAx(4));
semilogy(fsl(2,:), KfMAP./fsl, '-'); hold on;
semilogy(fsl(2,:),   Kf0./fsl, ':'); hold off;
ylabel('volume diffusivity [m2/Pa s]');
set(gca,'YAxisLocation','right');

axes(hAx(5));
semilogy(fsl(2,:), fsl.^2./CvMAP, '-'); hold on;
semilogy(fsl(2,:), fsl.^2./Cv0  , ':'); hold off;
xlabel('liquid fraction');
ylabel('segregation coefficient [m2/Pa s]');

axes(hAx(6));
semilogy(fsl(2,:), fsl.^2./CfMAP, '-'); hold on;
semilogy(fsl(2,:), fsl.^2./Cf0  , ':'); hold off;
xlabel('liquid fraction');
ylabel('compaction coefficient [Pa s]');
set(gca,'YAxisLocation','right');

annotation('textbox','Position',[0.4,0.88,0.2,0.1],'FitBoxToText','on',...
    'String','solid: MAP model from CATMIP, dotted: previous calib',...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',16);
SaveFigure([rundir runID '_comparexMAPtopreviouscalib.pdf']);

%% save outputs

save([rundir runID '_par.mat'],'PHS','runID','distrib','Niter','Nstep','pllopt','vname',...
    'rho0','eta0','d0','xccep','xd','ftot','data','dcat','sigm','xout','Pout');
diary off


