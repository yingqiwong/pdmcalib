function [xregm, regmname] = flowregime_twophase (fl, dsc, d0)

dscsl = squeeze(dsc(1,2,:));

[dscmax,poroind] =  max(dscsl)
        suspind  = find(dscsl > 10*d0(2), 1 );

regmname = {'POROUS','MUSH','SUSPENSION'};
xregm = [0, fl(poroind); fl(poroind), fl(suspind); fl(suspind), 1.0000];


end