%
% script to plot various functions to calibrate three-phase model
% YQW, 23 March 2022

clear all; clc; %close all; 
Addpaths

fprintf('\n\nThe list of available calibration functions are: \n');
dir('../src/calibfuncs/*.m')
fprintf('\n\n');

load ocean.mat

%% phase properties

PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+4;1e-3];   % pure-phase viscosities
d0   = [1e-3 ;1e-3;1e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
% original
A = [   0.60, 0.25, 0.30;
        0.20, 0.20, 0.20;
        0.20, 0.20, 0.20; ];  % permission slopes
B = [   0.30, 0.15, 0.55;
        0.48, 0.02, 0.50;
        0.80, 0.08, 0.12; ];  % permission step locations
C = [   0.20, 0.20, 0.20;
        0.60, 0.60, 0.12;
        0.20, 0.25, 0.50; ];  % permission step widths

% for fiddling
A = [   0.60, 0.25, 0.30;
        0.20, 0.20, 0.20;
        0.20, 0.20, 0.20; ];  % permission slopes
B = [   0.30, 0.15, 0.55;
        0.48, 0.02, 0.50;
        0.47, 0.48, 0.05; ];  % permission step locations
C = [   0.20, 0.20, 0.20;
        0.60, 0.60, 0.12;
        0.07, 0.02, 0.20; ];  % permission step widths

Nf = 200;       % number of points to generate data on

%% plot coefficients

f = SetUp3PhsMatrix(Nf);            % three-phase variation
[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
Plot3PhaseCoeff(f, cat(3,f.^2./Cv,f.^2./Cf), 'scl', 'log', 'PHS', PHS, 'cfname', {'(\phi^2/C_v)','(\phi^2/C_\phi)'});

%% parmigiani critical gas fraction and permeability

f = SetUp3PhsMatrix(Nf);            % three-phase variation
[~, ~, ~, ~, ~, Xf] = SegCompLength(f, eta0, d0, A, B, C);

[kmvp, fpm] = gcrit_parmi(f);
[~, ~, ~, Cvpm] = SegCompLength(fpm, eta0, d0, A, B, C);
kmodel = fpm(3,:).^2./Cvpm(3,:).*eta0(3);

% plot
figure('Name','parmigiani model'); 
hAx = setupaxes(1,3,'left',2,'gapw',2.3,'right',2);
fontsize = get(gcf,'defaultaxesfontsize');

axes(hAx(1)); 
plot(fpm(1,:), fpm(3,:)); grid on;
xlabel('crystal fraction'); ylabel('critical mvp fraction');

axes(hAx(2)); 
semilogy(fpm(1,:), kmvp, fpm(1,:), kmodel); grid on
xlabel('crystal fraction'); ylabel('mvp permeability [m$^2$]');
legend('parmigiani','closure model','Location','best');

axes(hAx(3)); 
hAx(3).Position(3) = 1.2*hAx(3).Position(3);
ternsurf(f(1,:), f(2,:), f(3,:), squeeze(Xf(3,3,:)));
colormap(ocean);
hold on;
hpm = ternplot(fpm(1,:),fpm(2,:),fpm(3,:), 'w', 'linewidth', 4);
hpm.ZData = 2*ones(size(hpm.XData));
shading interp; view([0,90]);
vertexlabel(PHS{1},PHS{2},PHS{3});
cb = colorbar('eastoutside','TickLabelInterpreter','latex');
cb.Position(1) = cb.Position(1) + 0.002;
hAx(3).Position(1) = hAx(3).Position(1) - 0.02;
title(cb, '$X_\phi^{mvp-mvp}$','Interpreter','latex');

sgtitle('Outgassing, Parmigiani 2017 model', 'FontSize', fontsize);

% SaveFigure('Figures/parmigiani_gcrit_trialerrorfit');

%% three-phase viscosity models

f = SetUp3PhsMatrix(Nf);            % three-phase variation
[~, Kv] = SegCompLength(f, eta0, d0, A, B, C);

eta_bb    = visc_birnbaum(f, eta0(2));
eta_truby = visc_truby(   f, eta0(2));

% plot
figure('Name','three-phase viscosities'); 
hAx = setupaxes(1,3,'height',9,'top',1.8,'bot',0.5,'gapw',0.2,'left',0.2,'right',0.2);
fontsize = get(gcf,'defaultaxesfontsize');

axes(hAx(1));
ternpcolor(f(1,:), f(2,:), f(3,:), eta_bb, 10);
colormap(ocean);
vertexlabel('$\phi^1$','$\phi^2$','$\phi^3$'); 
set(gca,'ColorScale','log'); clims = caxis; 
cb = colorbar('Location','southoutside'); cb.TickLabelInterpreter = 'latex';
text(0.5,1,'Mixture viscosity, Birnbaum 2021','HorizontalAlignment','center','FontSize',fontsize);

axes(hAx(2));
ternpcolor(f(1,:), f(2,:), f(3,:), eta_truby, 10);
colormap(ocean);
vertexlabel('$\phi^1$','$\phi^2$','$\phi^3$'); 
set(gca,'ColorScale','log'); caxis(clims); 
cb = colorbar('Location','southoutside');  cb.TickLabelInterpreter = 'latex'; 
text(0.5,1,'Mixture viscosity, Truby 2015','HorizontalAlignment','center','FontSize',fontsize);

axes(hAx(3));
ternpcolor(f(1,:), f(2,:), f(3,:), sum(Kv,1), 10);
colormap(ocean);
vertexlabel('$\phi^1$','$\phi^2$','$\phi^3$'); 
set(gca,'ColorScale','log'); caxis(clims); 
cb = colorbar('Location','southoutside');  cb.TickLabelInterpreter = 'latex'; 
text(0.5,1,'Mixture viscosity, model','HorizontalAlignment','center','FontSize',fontsize);

% SaveFigure('Figures/coeffs_mixtureviscosity');

%% models on the solid-liquid axis (reduce to two-phase case)

fsl = SetUp3PhsMatrix(Nf, 3, 0.001);  
[~, Kv, Kf, Cv, Cf, ~] = SegCompLength(fsl, eta0, d0, A, B, C);

% viscosity models
[eta_er         ] = visc_er   (fsl, eta0(2));
[eta_co, zeta_co] = visc_costa(fsl, eta0(2));
[eta_ru, zeta_ru] = visc_rudge(fsl, eta0(1));
[eta_hk         ] = visc_hk   (fsl, eta0(1));
[eta_th, zeta_th] = visc_th   (fsl, eta0(1));

% volume diffusivity
[kphi_hs, Cseg_hs] = voldiff_hs   (fsl(2,:), fsl(1,:), eta0(2), d0(1));
[kphi_sg         ] = voldiff_segre(fsl(2,:), fsl(1,:), eta0(2), d0(2));
[         Cseg_kc] = perm_kc      (fsl(1,:), fsl(2,:), eta0(2), d0(1));

% plot
figure('Name','viscosities on solid-liquid axis'); 
hAx = setupaxes(2,2,'left',2,'gapw',2.5,'height',6.5);
set(gcf,'defaultaxescolororder', [[0,0,0];[0,0,0];lines(7)]); 

axes(hAx(1));
semilogy(fsl(2,:), Kv(1,:)./fsl(1,:), '-'); hold on;
semilogy(fsl(2,:), Kv(2,:)./fsl(2,:), ':'); 
semilogy(fsl(2,:), eta_co, '--'); 
semilogy(fsl(2,:), eta_ru); 
semilogy(fsl(2,:), eta_th); 
semilogy(fsl(2,:), eta_hk, '--'); 
semilogy(fsl(2,:), eta_er); hold off; 
ylim([1e0,1e20]);
legend('model, solid','model, liquid','Costa 2009','Rudge 2018','Takei \& Holtzman 2009',...
    'Hirth \& Kohlstedt 2003','Einstein-Roscoe','box','off');
xlabel('liquid fraction');
ylabel('mixture viscosity [Pa s]');

axes(hAx(2));
semilogy(fsl(2,:), Kf(1,:)./fsl(1,:), '-'); hold on;
semilogy(fsl(2,:), Kf(2,:)./fsl(2,:), ':');
semilogy(fsl(2,:), kphi_hs); 
semilogy(fsl(2,:), kphi_sg); hold off;
xlim([0,1]); ylim([1e-25,1e-5]);
legend('model, solid','model, liquid','hindered-Stokes','Segre 2001','Location','southeast','box','off');
xlabel('liquid fraction');
ylabel('volume diffusivity [m$^2$/Pa s]');

axes(hAx(3));
semilogy(fsl(2,:), fsl(1,:).^2./Cv(1,:), '-'); hold on;
semilogy(fsl(2,:), fsl(2,:).^2./Cv(2,:), ':'); 
semilogy(fsl(2,:), Cseg_hs); 
semilogy(fsl(2,:), Cseg_kc); hold off;
ylim([1e-25,1e-5]);
legend('model, solid','model, liquid','hindered-Stokes','Kozeny-Carman','Location','south','box','off');
xlabel('liquid fraction');
ylabel('segregation coefficient [m$^2$/Pa s]');

axes(hAx(4));
semilogy(fsl(2,:), fsl(1,:).^2./Cf(1,:), '-'); hold on;
semilogy(fsl(2,:), fsl(2,:).^2./Cf(2,:), ':'); 
semilogy(fsl(2,:), zeta_co, '--'); 
semilogy(fsl(2,:), zeta_ru);
semilogy(fsl(2,:), zeta_th); hold off;
ylim([1e0,1e20]);
legend('model, solid','model, liquid','Costa 2009','Rudge 2018','Takei \& Holtzman 2009','box','off');
xlabel('liquid fraction');
ylabel('compaction viscosity [Pa s]');

sgtitle('coefficients on solid-liquid axis');

% SaveFigure('Figures/coeffs_solidliquidaxis');

%% models on gas-solid axis

fgs = SetUp3PhsMatrix(Nf, 2, 0.001);  
[~, Kv, Kf, Cv, Cf] = SegCompLength(fgs, eta0, d0, A, B, C);

segKC = perm_kc(fgs(1,:), fgs(3,:), eta0(3), d0(1));
etahf = visc_haff(fgs);

% plot
figure('Name','coefficients on gas-solid axis'); 
hAx = setupaxes(2,2,'left',2,'gapw',2.5,'height',6.5);
set(gcf,'defaultaxescolororder', [[0,0,0];[0,0,0];lines(7)]); 

axes(hAx(1));
semilogy(fgs(3,:), Kv(3,:)./fgs(3,:), '-'); hold on;
semilogy(fgs(3,:), Kv(1,:)./fgs(1,:), ':'); axis manual
semilogy(fgs(3,:), eta0(3).*etahf, '--'); 
hold off; 
legend('model, gas','model, solid','haff model','box','off','Location','southwest');
xlabel('gas fraction');
ylabel('mixture viscosity [Pa s]');

axes(hAx(2));
semilogy(fgs(3,:), Kf(3,:)./fgs(3,:), '-'); hold on;
semilogy(fgs(3,:), Kf(1,:)./fgs(1,:), ':'); 
xlim([0,1]); 
legend('model, gas','model, solid','hindered-Stokes','Location','south','box','off');
xlabel('gas fraction');
ylabel('volume diffusivity [m$^2$/Pa s]');

axes(hAx(3));
semilogy(fgs(3,:), fgs(3,:).^2./Cv(3,:), '-'); hold on;
semilogy(fgs(3,:), fgs(1,:).^2./Cv(1,:), ':');
semilogy(fgs(3,:), segKC, '--'); hold off;
legend('model, gas','model, solid','Kozeny-Carman','Location','south','box','off');
xlabel('gas fraction');
ylabel('segregation coefficient [m$^2$/Pa s]');

axes(hAx(4));
semilogy(fgs(3,:), fgs(3,:).^2./Cf(3,:), '-'); hold on;
semilogy(fgs(3,:), fgs(1,:).^2./Cf(1,:), ':'); hold off;
legend('model, gas','model, solid','box','off');
xlabel('gas fraction');
ylabel('compaction viscosity [Pa s]');

SaveFigure('Figures/coeffs_gassolidaxis');


%% models on gas-liquid axis

fgl = SetUp3PhsMatrix(Nf, 1, 0.001);  
[~, Kv, Kf, Cv, Cf] = SegCompLength(fgl, eta0, d0, A, B, C);

% various bubble-liquid viscosity models at high capillary number
eta_dr =    visc_dr(fgl, eta0(2));
eta_bd =    visc_bd(fgl, eta0(2));
eta_ll =    visc_ll(fgl, eta0(2));
eta_md = visc_mader(fgl, eta0(2), 10);

% hindered Stokes on dispersed bubble phase
[kphi_hs,   seg_hs] = voldiff_hs(fgl(2,:), fgl(3,:), eta0(2), d0(3));

% plot
figure('Name','coefficients on gas-liquid axis'); 
hAx = setupaxes(2,2,'left',2,'gapw',2.5,'height',6.5);
set(gcf,'defaultaxescolororder', [[0,0,0];[0,0,0];lines(7)]); 

axes(hAx(1));
semilogy(fgl(3,:), Kv(3,:)./fgl(3,:), '-'); hold on;
semilogy(fgl(3,:), Kv(2,:)./fgl(2,:), ':'); axis manual
semilogy(fgl(3,:), eta_dr, '--');
semilogy(fgl(3,:), eta_bd, '--');
semilogy(fgl(3,:), eta_ll, '--');
semilogy(fgl(3,:), eta_md, '--');
hold off; 
legend('model, gas','model, liquid','Ducamp Raj','Bagdassarov Dingwell',...
    'Llewellin','Mader, Ca=10','box','off','Location','southwest');
xlabel('gas fraction');
ylabel('mixture viscosity [Pa s]');

axes(hAx(2));
semilogy(fgl(3,:), Kf(3,:)./fgl(3,:), '-'); hold on;
semilogy(fgl(3,:), Kf(2,:)./fgl(2,:), ':'); 
semilogy(fgl(3,:), kphi_hs, '--');  hold off;
xlim([0,1]); 
legend('model, gas','model, liquid','hindered-Stokes','Location','south','box','off');
xlabel('gas fraction');
ylabel('volume diffusivity [m$^2$/Pa s]');

axes(hAx(3));
semilogy(fgl(3,:), fgl(3,:).^2./Cv(3,:), '-'); hold on;
semilogy(fgl(3,:), fgl(2,:).^2./Cv(2,:), ':');
semilogy(fgl(3,:), seg_hs, '--');  hold off;
legend('model, gas','model, liquid','hindered-Stokes','Location','south','box','off');
xlabel('gas fraction');
ylabel('segregation coefficient [m$^2$/Pa s]');

axes(hAx(4));
semilogy(fgl(3,:), fgl(3,:).^2./Cf(3,:), '-'); hold on;
semilogy(fgl(3,:), fgl(2,:).^2./Cf(2,:), ':'); hold off;
legend('model, gas','model, liquid','box','off');
xlabel('gas fraction');
ylabel('compaction viscosity [Pa s]');

sgtitle('coefficients on gas-liquid axis');

% SaveFigure('Figures/coeffs_gasliquidaxis');
