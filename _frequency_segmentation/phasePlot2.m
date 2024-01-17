    function h =  phasePlot2(curv, K, dKdt, s, ax, xl, yl)
    
        Kc_slice   = interp1(curv, K',   s, 'makima')';
        dKdc_slice = interp1(curv, dKdt',s, 'makima')';
        h = plot(ax, Kc_slice, dKdc_slice, 'LineWidth', 2);
        xm = abs(min(xl));
        ym = abs(min(yl));
        xlim([-xm xm])
        ylim([-ym ym])
        xlabel('K')
        ylabel('dKdt')
        title(sprintf('Current position: %2.0f%%', s))
        set(gca, 'FontSize', 12)
    end
