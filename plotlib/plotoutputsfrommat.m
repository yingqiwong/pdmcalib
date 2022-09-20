function [] = plotoutputsfrommat (folder, RunID, save_opt)

load([folder, RunID, '/', RunID, '_par.mat']);

% posterior distributions
[xMAP] = plotdistribs(xout, Pout, vname, xd, distrib);

if save_opt, SaveFigure([rundir runID '_posteriordistributions.pdf']); end

% corner plots to see correlations between posterior distributions
figs = plotcorrelations(xout, Pout, vname, xd, distrib);
if save_opt
for fi = 1:length(figs)
    SaveFigure([rundir runID '_correlations_' num2str(fi) '.pdf'], figs(fi));
end
end

% print out MAP model
[A_MAP, B_MAP, C_MAP] = permvec2mat(xMAP, 'log')

% plot the fit to the datasets
plotfittodata(ftot, Ns, data, dhf(xMAP), dcat, sigm);
if save_opt, SaveFigure([rundir runID '_xMAP_fittodata.pdf']); end

% now calculate connectivity and coefficients for plotting
fout = SetUp3PhsMatrix(200);
[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(fout, eta0, d0, A_MAP, B_MAP, C_MAP);

% plot permissions
ax = Plot3PhasePerm(fout, Xf, PHS);
if save_opt, SaveFigure([rundir runID '_xMAP_connectivity.pdf']); end

% plot effective viscosity and volume diffusivity
Plot3PhaseCoeff(fout,cat(3,Kv,Kf),'scl','log','PHS',PHS,'cfname',{'(K_v/\phi)','(K_\phi/\phi)'});
if save_opt, SaveFigure([rundir runID '_xMAP_fluxcoeff.pdf']); end

% plot segregation and compaction coefficients
Plot3PhaseCoeff(fout,cat(3,fout.^2./Cv,fout.^2./Cf),'scl','log','PHS',PHS,'cfname',{'(\phi^2/C_v)','(\phi^2/C_\phi)'});
if save_opt, SaveFigure([rundir runID '_xMAP_segcompcoeff.pdf']); end

% plot weights
omCv = Cv./sum(Cv,1);
omCf = Cf./sum(Cf,1);
Plot3PhaseCoeff(fout,cat(3,omCv,omCf),'PHS',PHS,'cfname',{'(\omega_{Cv})','(\omega_{C\phi})'});
if save_opt, SaveFigure([rundir runID '_xMAP_vpweights.pdf']); end


end