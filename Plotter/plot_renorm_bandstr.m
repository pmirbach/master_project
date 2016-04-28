function [ren_bandstr_surf, ren_bandstr_path] = plot_renorm_bandstr(Ctrl,Para,k,Daten,layout)

if any(Ctrl.plot.ren_bs)
    Daten_plot = Daten * 1e-3;
    titlestr = {'V \uparrow','L1 \uparrow','L2 \uparrow',...
        'V \downarrow','L1 \downarrow','L2 \downarrow'};
    titlestr = [titlestr, titlestr];
end

if Ctrl.plot.ren_bs(1) == 1
    if Ctrl.plot.entireBZ == 1
        k_plot = k;
    else
        k_plot = k(:,:,1);
    end
    ren_bandstr_surf = plot_surf(Para,k_plot,Daten_plot,layout);
    for ii = 1:size(ren_bandstr_surf,2)
        figure(ren_bandstr_surf)
        ax = flip(findobj(ren_bandstr_surf(ii),'type','axes'),1);
        cb = flip(findobj(ren_bandstr_surf(ii),'type','colorbar'),1);
        
        for jj = 1:size(ax,1)
            ax(jj).FontSize = 12;
            ax(jj).Title.String = titlestr{jj};
            ax(jj).Position = ax(jj).Position  + [-.03 -.03 0.03 0.03];
            cb(jj).FontSize = 10;
            cb(jj).Title.String = 'E - E_f (eV)';
        end
        mtit('Dispersion der Energiebänder','fontsize',18,...
            'xoff',0,'yoff',.03)
    end
    if Ctrl.plot.save == 1
        for ii = 1:size(ren_bandstr_surf,2)
            filename = sprintf(['/home/pmirbach/Masterarbeit/MoS2/',...
                'Anlauf2/Ergebnisse/bands_2d_sub%d'], ii);
            export_fig(h_band_2d(ii),filename, '-png');
        end
    end
else
    ren_bandstr_surf = [];
end

if Ctrl.plot.ren_bs(2) == 1
    ren_bandstr_path = plot_path(Ctrl,Para,k,Daten_plot,200);
       
    title_string = sprintf('Renormierte Bandstruktur bei D_0 = %g', Ctrl.carrier_density);
    
    ax = findobj(ren_bandstr_path,'type','axes');
    ax.Title.String = title_string;
    legend(titlestr,'location','northeastoutside')
    ax.YLabel.String = 'E - E_f (eV)';
    
    plot_lines = findobj(ren_bandstr_path,'Type','line');
    plot_lines_del = zeros(size(plot_lines));
    
    for ii = 1:size(plot_lines,1)
        if isequal(plot_lines(ii).Color,[0.5 0.5 0.5])
            plot_lines_del(ii) = 1;
        end
    end
    plot_lines(logical(plot_lines_del)) = [];
    
    
    for ii = 1:size(plot_lines,1)
        if ii <= size(plot_lines,1)/2
            plot_lines(ii).LineStyle = '-';
        else
            plot_lines(ii).LineStyle = ':';
            plot_lines(ii).LineWidth = 1;
        end
    end
    
    
    if Ctrl.plot.save == 1
        export_fig(['/home/pmirbach/Masterarbeit/MoS2/Anlauf2/', ...
            'Ergebnisse/bands_path'], '-png');
    end
else
    ren_bandstr_path = [];
end




