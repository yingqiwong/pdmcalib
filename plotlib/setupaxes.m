
function [hAx,pos] = setupaxes (Nrow, Ncol, varargin)
% 
% [hAx,pos] = setupaxes (Nrow, Ncol, varargin)
% 
% function to set up axis with specific commands

ax.height = 6.00; 
ax.bot    = 2.00; 
ax.top    = 2.00;
ax.gaph   = 2.00; 

ax.width  = 9.00; % paper width
ax.left   = 2.00;
ax.right  = 2.00;
ax.gapw   = 2.00;

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    ax.(args{1,ia}) = args{2,ia};
end

fh = ax.bot  + Nrow*ax.height + (Nrow-1)*ax.gaph + ax.top;
fw = ax.left + Ncol*ax.width  + (Ncol-1)*ax.gapw + ax.right; 

set(gcf,'Units','centimeters','Position',[5 5 fw fh]);
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(gcf,'defaultaxesfontsize',14);
set(gcf,'defaultlinelinewidth',1);
set(gcf,'Resize','off');

if Nrow>1 || Ncol>1
    % predefine axis positions
    pos = zeros(Nrow*Ncol, 4);
    
    for ri = 1:Nrow
        for ci = 1:Ncol
            xc = ax.left + (   ci-1)*ax.width  +   (ci-1)*ax.gapw;
            yc = ax.bot  + (Nrow-ri)*ax.height +(Nrow-ri)*ax.gaph; 
             
            pos((ri-1)*Ncol + ci,:) = [xc, yc, ax.width, ax.height];
            hAx((ri-1)*Ncol + ci  ) = axes('Units','Centimeters','position',pos((ri-1)*Ncol + ci,:) );
        end
    end
    
elseif Nrow==1 && Ncol==1
    
    hAx = axes('Units','centimeters','Position',[ax.left,ax.bot,ax.width,ax.height]);
    
end


end
