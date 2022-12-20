

%%  calculate connectivity and coefficients

[dsc , KvMAP, KfMAP, CvMAP, CfMAP, XfMAP] = SegCompLength(fsl, eta0, d0, A_MAP, B_MAP, C_MAP);
[dsc0, Kv0  , Kf0  , Cv0  , Cf0  , Xf0  ] = SegCompLength(fsl, eta0, d0, A0   , B0   , C0   );

%% plot connectivity and coefficients

figure;
hAx = setupaxes(4,2,'gaph',0.8,'gapw',0.3);
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
ylabel('segregation coefficient [m2/Pa s]');

axes(hAx(6));
semilogy(fsl(2,:), fsl.^2./CfMAP, '-'); hold on;
semilogy(fsl(2,:), fsl.^2./Cf0  , ':'); hold off;
ylabel('compaction coefficient [Pa s]');
set(gca,'YAxisLocation','right');

axes(hAx(7));
plot(fsl(2,:), CvMAP./sum(CvMAP,1), '-'); hold on;
plot(fsl(2,:), Cv0  ./sum(Cv0  ,1), ':'); hold off;
xlabel('liquid fraction');
ylabel('velocity weights');

axes(hAx(8));
plot(fsl(2,:), CfMAP./sum(CfMAP,1), '-'); hold on;
plot(fsl(2,:), Cf0  ./sum(Cf0  ,1), ':'); hold off;
xlabel('liquid fraction');
ylabel('pressure weights]');
set(gca,'YAxisLocation','right');

annotation('textbox','Position',[0.4,0.9,0.2,0.1],'FitBoxToText','on',...
    'String','solid: MAP model from CATMIP, dotted: previous calib',...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',16,'Interpreter','latex');
SaveFigure([rundir runID '_comparexMAPtopreviouscalib.pdf']);

%% plot segregation-compaction length

[~,poroind] =  max(squeeze( dsc(1,2,:) ) );
   suspind  = find(squeeze( dsc(1,2,:) ) > 10*d0(2), 1 );

regmname = {'POROUS','MUSH','SUSPENSION'};
xregm = [0, fsl(2,poroind); fsl(2,poroind), fsl(2,suspind); fsl(2,suspind), 1.0000]
tmpcolors  = lines(3);
rectcolors = [[tmpcolors(1,:); tmpcolors(3,:); tmpcolors(2,:)], 0.12*ones(3,1)];

figure;
setupaxes(1,1,'top',0.8,'bottom',-1,'right',0.5,'height',12);
hsc(1) = semilogy(fsl(2,:), squeeze(dsc(1,2,:)), 'k-');
hold on; 
hsc(3) = semilogy(fsl(2,:), squeeze(dsc0(1,2,:)), 'k:'); 
hsc(2) = semilogy(fsl(2,:), squeeze(dsc (2,1,:)), '-', 'Color', 0.7*ones(1,3));
hsc(4) = semilogy(fsl(2,:), squeeze(dsc0(2,1,:)), ':', 'Color', 0.7*ones(1,3)); 
plot(xlim, d0(1)*ones(1,2), 'k--', 'linewidth', 1); 
hold off;
text(0.1, d0(1), 'grain size', 'FontSize', 12, 'VerticalAlignment', 'bottom');

% plot regimes
ylimits = [1e-7,1e3]; ylim(ylimits);
for ri = 1:3
    r = rectangle('Position',[xregm(ri,1),ylimits(1),diff(xregm(ri,:)),diff(ylimits)], ...
        'FaceColor',rectcolors(ri,:),'EdgeColor','none');
    uistack(r,'bottom');
    text(xregm(ri,1)+0.5*diff(xregm(ri,:)),1,regmname{ri},...
        'FontSize',12,'Rotation',90,'HorizontalAlignment','center','Color',rectcolors(ri,1:3));
end

set(gca,'XTick', 0:0.2:1);
set(gca,'YTick', 10.^(-10:2:6),'YMinorTick','off');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('$\delta_{0}^{jk}$ [m]');
title('segregation-compaction length');
legend(hsc,...
    '$\delta_{0}^{\ell s}$, this inversion', '$\delta_{0}^{s \ell}$, this inversion', ...
    '$\delta_{0}^{\ell s}$, initial model' , '$\delta_{0}^{s \ell}$, initial model',...
    'Location', 'southoutside','NumColumns',2,'box','off');

SaveFigure([rundir runID '_segcomplength.pdf']);

%% plot velocity and pressure scales

df   = 0.01;
g    = 9.81;
D    = [1; 10; 100];
drhobar = abs(diff(rho0));
drhoseg = abs(rho0 - sum(fsl.*rho0,1));

uSegr =    drhoseg.*g.*fsl.^2./CvMAP;
uBar  = df*drhobar.*g.*D.^2  ./sum(KvMAP,1);
pComp =    drhoseg.*g.*[squeeze(dsc(1,2,:))'; squeeze(dsc(2,1,:))'];
pBar  = df*drhobar.*g.*D.*ones(1,Nf);

figure; 
hAx = setupaxes(1,2,'top',0.8,'right',0.5,'height',8);

axes(hAx(1)); 
semilogy(fsl(2,:), uSegr); hold on;
semilogy(fsl(2,:), uBar, 'Color', 0.7*ones(1,3)); hold off;
text(0.9*ones(3,1),uBar(:,1), num2str(D(:)), 'VerticalAlignment','top');
ylimits = ylim;
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('Velocity scale [m/s]');

for ri = 1:3
    r = rectangle('Position',[xregm(ri,1),ylimits(1),diff(xregm(ri,:)),diff(ylimits)], ...
        'FaceColor',rectcolors(ri,:),'EdgeColor','none');
    uistack(r,'bottom');
end


axes(hAx(2)); 
semilogy(fsl(2,:), pComp); hold on;
semilogy(fsl(2,:), pBar, 'Color', 0.7*ones(1,3)); hold off;
text(0.1*ones(3,1),pBar(:,end), num2str(D(:)));
ylimits = ylim;
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('Pressure scale [Pa]');

for ri = 1:3
    r = rectangle('Position',[xregm(ri,1),ylimits(1),diff(xregm(ri,:)),diff(ylimits)], ...
        'FaceColor',rectcolors(ri,:),'EdgeColor','none');
    uistack(r,'bottom');
end

SaveFigure([rundir runID '_vpscales.pdf']);
