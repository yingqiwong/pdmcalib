
function [xMAP] = plotdistribs (x, p, vname, xd, distrib)
% 
% [xMAP] = plotdistribs (x, p, vname, xbnds, xd, distrib)
% 
% plots output distributions from the output of catmip or mcmc
% 
% INPUTS
% x         matrix of posterior-sampled model parameters [Niter x Nvar]
% p         output posterior probabilities [Niter x 1]
% vname     model parameter names [Nvar x 1]
% distrib   name of distribution: 'uniform' or 'normal'
% 
% OUTPUTS
% xMAP      maximum a posteriori model [1 x Nvar]
% 
% YQW, 20 April 2022

% collect info on some sizes
Nbins = 100;
NPHS  = sqrt(length(vname)/3);
Nx    = size(x,1);

% get maximum a posteriori model
[~, xMAPind] = max(p);
xMAP = x(xMAPind,:);

% get indices of x corresponding to A, B, C parameters
Ai = find(contains(vname, 'A'));
Bi = find(contains(vname, 'B'));
Ci = find(contains(vname, 'C'));

switch distrib
    case 'normal',  xbnds = xd(:,1) + 4*[-1,1].*xd(:,2);    % limits are 4 standard deviations
    case 'uniform', xbnds = xd;                             % limits are parameter bounds
end


figure;
hAx = setupaxes(3,NPHS^2,'width',5,'height',4,'gapw',1,'left',0.5,'right',0.4,'top',1,'bot',1.3);
axi = 1;
for mi = 1:NPHS
    for ni = 1:NPHS
        pind = mi + (ni-1)*NPHS;
        
        axes(hAx(0*NPHS^2+axi));
        ha = histogram(x(:,Ai(pind)),Nbins,'EdgeColor','none'); hold on;
        plot(xMAP(Ai(pind)  )*ones(1,2), ylim, 'r:');           hold off;
        xlim(xbnds(Ai(pind),:));
        set(gca,'YTickLabel',[],'XTickLabel',num2str(10.^get(gca,'XTick')'));
        title(vname{Ai(pind)});
        
        axes(hAx(1*NPHS^2+axi));
        hb = histogram(x(:,Bi(pind)),Nbins,'EdgeColor','none'); hold on;
        plot(xMAP(Bi(pind)  )*ones(1,2), ylim, 'r:');           hold off;
        xlim(xbnds(Bi(pind),:));
        set(gca,'YTickLabel',[],'XTickLabel',num2str(10.^get(gca,'XTick')'));
        title(vname{Bi(pind)});
        
        axes(hAx(2*NPHS^2+axi));
        hc = histogram(x(:,Ci(pind)),Nbins,'EdgeColor','none'); hold on;
        plot(xMAP(Ci(pind)  )*ones(1,2), ylim, 'r:');           hold off;
        xlim(xbnds(Ci(pind),:));
        set(gca,'YTickLabel',[],'XTickLabel',num2str(10.^get(gca,'XTick')'));
        title(vname{Ci(pind)});
        
        % plot normal distribution prior
        if strcmp(distrib, 'normal')
            axes(hAx(0*NPHS^2+axi)); plotnormd(xd(Ai(pind),1), xd(Ai(pind),2), Nx*ha.BinWidth);
            axes(hAx(1*NPHS^2+axi)); plotnormd(xd(Bi(pind),1), xd(Bi(pind),2), Nx*hb.BinWidth);
            axes(hAx(2*NPHS^2+axi)); plotnormd(xd(Ci(pind),1), xd(Ci(pind),2), Nx*hc.BinWidth);
        end
        
        axi = axi + 1;
        
    end
end

end

function [] = plotnormd (mu, sigma, N)
% plot normal distribution prior

xlimits = xlim;
xvec    = linspace(xlimits(1), xlimits(2), 1001);

hold on;
plot(xvec, sigma*sqrt(2*pi)*N*normpdf(xvec, mu, sigma), 'r-');
hold off;

end
