function subdivisionDiagram(f, levels)
% SUBDIVISIONDIAGRAM draw a diagram to show subdivision and polynomial
% degrees.
g = chebfun(f);
LW = 'linewidth'; lw = 2;
for lvl = 0:levels
    plot([-1 1],[lvl lvl],'k',LW,lw), hold on,
    xx = linspace(-1,1,2^lvl+1);
    newg = chebfun(f, xx); j = 1;
    for x = xx(1:end-1)
        plot([x x],[-.1 .1]+lvl,'k',LW,lw)
        funs = newg.funs;
        text(x+2.^(-lvl)-.09,lvl+.1,sprintf('N = %u',length(funs{j})-1));
        j = j + 1;
    end
    plot([xx(end) xx(end)],[-.1 .1]+lvl,'k',LW,lw)
end
set(gca,'ytick',[])
set(gca,'yticklabel',[])
hold off
end
