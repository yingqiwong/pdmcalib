
function [] = plotfittodata (f, data, dhat, dcat, sigma)

dcatNames = unique(dcat,'stable');
Ndcat = length(dcatNames);
[Nrow, Ncol] = GetSubplotRowCol(Ndcat);

figure;
for dci = 1:Ndcat
    subplot(Nrow,Ncol,dci)
    fi = strcmp(dcat, dcatNames{dci});
    errorbar(f(2,fi), data(fi), sigma(fi), '.', 'MarkerSize', 20); 
    hold on;
    plot(f(2,strcmp(dcat,dcatNames{dci})), dhat(strcmp(dcat,dcatNames{dci})), '+'); 
    hold off;
    legend('observed','model-predicted','location','best');
    title(dcatNames{dci});
end

end

