

function [vout, fout] = samplevalidvalues (v, f, Ns, randflag)

f(:,isnan(v(1,:))) = [];
v(:,isnan(v(1,:))) = [];

if (randflag)
    % randomly sample among valid values
    ind = randperm(size(v,2), Ns);
else
    % sample evenly through valid values
    ind = unique(round(linspace(1,size(v,2),Ns)));
end

fout = f(:,ind);
vout = v(:,ind);

end

