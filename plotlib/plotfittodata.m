
function [] = plotfittodata (f, data, dhat, dcat, sigma)

dcatNames   = unique(dcat,'stable');
Ndcat       = length(dcatNames);

Nrow = floor(sqrt(Ndcat));
Ncol = ceil(Ndcat/Nrow);

figure;
hAx = setupaxes(Nrow, Ncol);

for di = 1:Ndcat
    axes(hAx(di));
    fi = strcmp(dcat, dcatNames{di});
    errorbar(f(2,fi), data(fi), sigma(fi), '.', 'MarkerSize', 20);
    hold on;
    plot(f(2,fi), dhat(fi), '+');
    hold off;
    legend('observed','model-predicted','location','best');
    title(dcatNames{di});
end

if length(hAx)>Ndcat
    for di = (Ndcat+1):length(hAx)
        hAx(di).Visible = 'off';
    end
end


end

