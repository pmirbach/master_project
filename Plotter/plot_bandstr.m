function [bandstr_surf, bandstr_path] = plot_bandstr(Ctrl,Parameter,k,Daten,layout)

if any(Ctrl.plot.tb)
    Daten_plot = Daten * 1e-3;
    titlestr = {'V \uparrow','L1 \uparrow','L2 \uparrow',...
        'V \downarrow','L1 \downarrow','L2 \downarrow'};
end

if Ctrl.plot.tb(1) == 1
    if Ctrl.plot.entireBZ == 1
        k_plot = k;
    else
        k_plot = k(:,:,1);
    end
    bandstr_surf = plot_surf(Parameter,k_plot,Daten_plot,layout);
    for ii = 1:size(bandstr_surf,2)
        figure(bandstr_surf)
        ax = flip(findobj(bandstr_surf(ii),'type','axes'),1);
        cb = flip(findobj(bandstr_surf(ii),'type','colorbar'),1);
        
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
        for ii = 1:size(bandstr_surf,2)
            filename = sprintf(['/home/pmirbach/Masterarbeit/MoS2/',...
                'Anlauf2/Ergebnisse/bands_2d_sub%d'], ii);
            export_fig(h_band_2d(ii),filename, '-png');
        end
    end
else
    bandstr_surf = [];
end

if Ctrl.plot.tb(2) == 1
    bandstr_path = plot_path(Ctrl,Parameter,k,Daten_plot,200);
    
    ax = findobj(bandstr_path,'type','axes');
    ax.Title.String = 'Bandstruktur entlang eines Pfades durch die BZ';
    legend(titlestr,'location','northeastoutside')
    ax.YLabel.String = 'E - E_f (eV)';
    
    if Ctrl.plot.save == 1
        export_fig(['/home/pmirbach/Masterarbeit/MoS2/Anlauf2/', ...
            'Ergebnisse/bands_path'], '-png');
    end
else
    bandstr_path = [];
end




