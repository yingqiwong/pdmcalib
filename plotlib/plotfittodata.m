
function [] = plotfittodata (f, Ns, data, dhat, dcat, sigma)

Ndcat   = floor(size(f,2)./Ns);
dcatvec = cell(Ndcat,1);

Nrow = floor(sqrt(Ndcat));
Ncol = ceil(Ndcat/Nrow);

figure;
hAx = setupaxes(Nrow, Ncol, 'gaph', 3);

for di = 1:Ndcat
    axes(hAx(di));
    
    fi = (di-1)*Ns + (1:Ns);
    dcatvec(di) = dcat(fi(1));
    
    % need to select which phase to plot on x axis (may change because
    % sometimes one phase is fixed
    switch dcatvec{di}
        case 'segr3'
            pi = 1;
        otherwise
            pi = 2;
    end
    
    errorbar(f(pi,fi), data(fi), sigma(fi), '.', 'MarkerSize', 20);
    hold on; plot(f(pi,fi), dhat(fi), '+');  hold off;
    
    xlabel(['phs ' num2str(pi)']);
    legend('observed','model-predicted','location','best');
    title(dcatvec{di});
end

if length(hAx)>Ndcat
    for di = (Ndcat+1):length(hAx)
        hAx(di).Visible = 'off';
    end
end


end

% possible options for dcat:
% etamix    mixture viscosity
% voldmix   mixture volume diffusion
% comp1     compaction coefficient of phase 1
% segr1     segregation coefficient of phase 1
% segr2     segregation coefficient of phase 2
% segr3     segregation coefficient of phase 3
% mvpcrt    critical mvp fraction for channelization onset
