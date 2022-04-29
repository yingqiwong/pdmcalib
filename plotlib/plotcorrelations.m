
function [fig] = plotcorrelations (x, p, vname, xd, distrib)
%
% [xMAP] = plotcorrelations (x, p, vname, xbnds, xd, distrib)
%
% plots correlations between distributions from the output of catmip or mcmc
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
NPHS = sqrt(length(vname)/3);
Nplt = 3*NPHS;
Nbin = 100;
Nx   = size(x,1);

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


% plot correlations
for iphs = 1:NPHS
    fig(iphs) = figure;
    hAx = setupaxes(NPHS^2,NPHS^2,'width',5,'height',4,'gapw',0,'gaph',0,'left',0.5,'right',1.5,'top',1,'bot',1.3);
    
    xind = (0:2)*NPHS^2 + iphs + (0:3:6)';
    
    vphs    = vname(  xind(:)   );
    xphs    =     x(:,xind(:)   );
    xMAPphs =  xMAP(  xind(:)   );
    xbndphs =  xbnds( xind(:), :);
    xdphs   =    xd(  xind(:), :);
    
    for mi = 1:Nplt
        for ni = 1:Nplt
            axes(  hAx((mi-1)*Nplt + ni) );
            
            if mi==ni
                % plot histogram
                hhist = histogram(xphs(:,mi), Nbin, 'EdgeColor', 'none');
                hold on;
                xlim(xbndphs(mi,:));

                % maximum a posteriori model
                plot(xMAPphs(mi)*ones(1,2), ylim, 'r:');
                hold off;
                
                % plot normal distribution?
                if strcmp(distrib, 'normal')
                    plotnormd(xdphs(mi,1), xdphs(mi,2), Nx*hhist.BinWidth);
                end
                
                xlabel(vphs{mi});
                
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',num2str(10.^get(gca,'XTick')'));
                
            elseif mi<ni
                dscatter(xphs(:,ni), xphs(:,mi),'plottype','contour');
                set(gca,'XTickLabel',[]);
                xlim(xbndphs(ni,:)); ylim(xbndphs(mi,:));
                if ni==Nplt
                    set(gca,'YAxisLocation','right');
                    set(gca,'YTickLabel',num2str(10.^get(gca,'YTick')','%.2f'));
                else
                    set(gca,'YTickLabel',[]);
                end
            else
                set(gca,'visible', 'off');
            end
            
        end
    end
    
    annotation('textbox','Position',[0.45,0.9,0.1,0.1],'String',...
        ['Correlation Plot for Phase ' num2str(iphs) ' (red line: MAP model)'],...
        'HorizontalAlignment','center','FontSize',20,'EdgeColor','none','Interpreter','latex');
end



end




