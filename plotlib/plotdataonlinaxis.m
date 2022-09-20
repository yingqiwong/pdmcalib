
function [] = plotdataonlinaxis (f, dcat)

dcatNames = unique(dcat,'stable');
Ndcat     = length(dcatNames);
colors    = lines(Ndcat);
mkr       = {'o','+','*','_','s','x','d','^','p','>','<'};

if Ndcat>length(mkr), mkr = repmat(mkr,1,20); end

figure;
setupaxes(1,1,'height',8,'width',12);
for dci = 1:Ndcat
    fi = strcmp(dcat, dcatNames{dci});
    htb(dci) = plot(f(1,fi),dci.*ones(length(f(1,fi)),1),...
        mkr{dci},'color',colors(dci,:),'markersize',5, 'linewidth',1);
    hold on;
end

ylim([0.5,Ndcat+0.5]);
set(gca,'YTick',1:Ndcat,'YTickLabel',dcatNames);
xlabel('solid fraction');

end