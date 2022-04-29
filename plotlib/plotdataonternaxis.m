
function [] = plotdataonternaxis (f, dcat)

dcatNames = unique(dcat,'stable');
Ndcat     = length(dcatNames);
colors    = lines(Ndcat);
mkr       = {'o','+','*','_','s','x','d','^','p','>','<'};

if Ndcat>length(mkr), mkr = repmat(mkr,1,20); end

figure;
setupaxes(1,1,'height',10,'width',15);
terngrid(10,0);
vertexlabel('$\phi^s$','$\phi^\ell$','$\phi^g$');
hold on;
for dci = 1:Ndcat
    fi = strcmp(dcat, dcatNames{dci});
    htb(dci) = ternplot(f(1,fi),f(2,fi),f(3,fi),mkr{dci},'color',colors(dci,:),'markersize',5, 'linewidth',1);
    htb(dci).ZData = 1.*ones(size(htb(dci).XData));
end
legend(htb,dcatNames(:),'location','eastoutside');

end