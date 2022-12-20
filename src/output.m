
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
plotfittodata(ftot, Ns, data, dhf(xMAP), dcat, sigm);
SaveFigure([rundir runID '_xMAP_fittodata.pdf']);

if NPHS == 2
    plot2phase;
elseif NPHS == 3
    plot3phase;
end

%% save outputs

save([rundir runID '_par.mat'],'PHS','runID','distrib','Niter','Nstep','pllopt','vname',...
    'rho0','eta0','d0','xccep','xd','ftot','data','dcat','sigm','xout','Pout');
diary off


