

function [hn] = plotnormd (mu, sigma, binwidth)
% overlay normal distribution prior on the histograms

xlimits = xlim;
xvec    = linspace(xlimits(1), xlimits(2), 1001);

hold on;
hn = plot(xvec, binwidth*normpdf(xvec, mu, sigma), 'r-');
hold off;

end


