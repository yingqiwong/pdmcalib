
% script to estimate some parameters in the costa 2009 model to fit
% material parameters
% YQW, 23 March 2022

clear all; clc; %close all; 
Addpaths

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
        0.48, 0.50, 0.02; ];  % permission step locations
C = [   0.20, 0.20, 0.20;
        0.60, 0.60, 0.12;
        0.07, 0.02, 0.20; ];  % permission step widths

%% initialise phase fractions

Nf  = 100;
fsl = SetUp3PhsMatrix(Nf, 3, 0.0001);   % along solid-liquid axis

%% first understand how parameters affect effective viscosity

[phistar, gamma, delta, xi] = olvbasparams;
params = [ phistar + [-0.2,0,0.2] ;
           gamma   + [-1  ,0,  1] ;
           delta   + [-10 ,0, 10] ;
           xi      * [ 0.1,1, 10] ];

figure;
hAx = setupaxes(2,2);

for pi = 1:size(params,2)
    axes(hAx(1));
    etatmp = costa_fit(fsl, eta0(2), params(1,pi), gamma, delta, xi);
    semilogy(fsl(2,:), etatmp); hold on;
    
    axes(hAx(2));
    etatmp = costa_fit(fsl, eta0(2), phistar, params(2,pi), delta, xi);
    semilogy(fsl(2,:), etatmp); hold on;
    
    axes(hAx(3));
    etatmp = costa_fit(fsl, eta0(2), phistar, gamma, params(3,pi), xi);
    semilogy(fsl(2,:), etatmp); hold on;
    
    axes(hAx(4));
    etatmp = costa_fit(fsl, eta0(2), phistar, gamma, delta, params(4,pi));
    semilogy(fsl(2,:), etatmp); hold on;
    
end

axes(hAx(1)); title('$\phi^*$'); legend(num2str(params(1,:)'));
axes(hAx(2)); title('$\gamma$'); legend(num2str(params(2,:)'));
axes(hAx(3)); title('$\delta$'); legend(num2str(params(3,:)')); xlabel('liquid fraction');
axes(hAx(4)); title('$\xi$'   ); legend(num2str(params(4,:)')); xlabel('liquid fraction');

SaveFigure('Figures/costa_varyparams');


%% fit costa parameter to new solid viscosity
% looks like we just need to adjust xi

eta_max = costa_fit(fsl(:,end-1), eta0(2), phistar, gamma, delta, xi);
xifit   = fsolve(@(xif) log10(costa_fit(fsl(:,end-1), eta0(2), phistar, gamma, delta, xif)) - log10(eta0(1)), xi);


%% plot fitted model against Einstein-Roscoe and Hirth-Kohlstsedt

eta    = costa_fit(fsl, eta0(2), phistar, gamma, delta, xifit);
eta_ER =   visc_er(fsl, eta0(2));
eta_HK =   visc_hk(fsl, eta0(1));
eta_tr =visc_truby(fsl, eta0(2));

figure;
setupaxes(1,1,'left',2,'width',10);
semilogy(fsl(2,:), eta_ER,  ':', 'linewidth', 4); hold on;
semilogy(fsl(2,:), eta_HK,  ':', 'linewidth', 4); 
semilogy(fsl(2,:), eta_tr,  ':', 'linewidth', 4);
semilogy(fsl(2,:), eta   , 'k-'); hold off;
legend('Einstein-Roscoe','Hirth-Kohlstedt','truby','fitted Costa','box','off');
xlabel('liquid fraction'); ylabel('viscosity [Pa s]');

SaveFigure('Figures/costa_fitted');


%% costa model

function eta =  costa_fit (f, eta_melt, phistar, gamma, delta, xi)

feff      = f(1:2,:)./sum(f(1:2,:),1);

% Costa+ effective mixture shear and compaction viscosities  (Costa et al., 2009)
B1        = 4.0;  % theoretical value = 2.5
phiscaled = feff(1,:)./phistar;
h         = (1-xi).*erf(sqrt(pi)./(2.*(1-xi)).*phiscaled.*(1+phiscaled.^gamma));

eta  = eta_melt .* (1+phiscaled.^delta) ./ (1-h).^(B1*phistar);

end

function [phistar, gamma, delta, xi] = olvbasparams ()

phistar   = 0.62;
gamma     = 3.25;
delta     = 24;
xi        = 4e-5;  % or 8e-8 (for B1 = 2.5)

end
    