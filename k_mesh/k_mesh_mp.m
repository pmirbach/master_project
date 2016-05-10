function [k, wk] = k_mesh_mp(Ctrl, Para)

qr = Ctrl.k_mesh_mp.qr;             % Feinheit des meshes

b1 = Para.k.GV(:,1);          % Reziproke Gittervektoren
b2 = Para.k.GV(:,2);


ur = (0:qr-1) / qr;                   % Einteilung des rez. Gittervektors
A = ones(numel(ur),1) * ur;          % Erzeugung eines meshes
k_mesh_x = A * b1(1) + A' * b2(1);
k_mesh_y = A * b1(2) + A' * b2(2);

k0 = [k_mesh_x(:), k_mesh_y(:)]';     % Alle k-Punkte


% Bestimmung der k-Punkte in der red BZ
k = pts_triangle_fun(k0, Para.BZred.symmpts{2}(:,1:3), 30 * eps);
% Bestimmung der Gewichte der k-Punkte
wk = pts_weight(k, Para.BZred.symmpts{2}(:,1:3), 30 * eps);

if round ( sum(wk) * Para.BZsmall.area * 1e3 ) ~= round( Para.BZ.area * 1e3 )
    warning off backtrace
    warning('Integrated weights do not agree with area of BZ!')
    warning on backtrace
end

% Erzeugung aller red. Dreiecke mit Spiegelungen und Drehungen
[k] = red_to_BZ(k);


%% Plots

if Ctrl.plot.k_mesh(1) == 1     % Plot über red. BZ mit den Gewichten
    figure
    hold on
    set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.6 0.9]);
    corners = [Para.BZred.symmpts{2} ...
        * [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], [0; 0]];
    
    in = wk == 6;
    on = wk == 3;
    symm = wk == 1;
    
    plot(corners(1,:),corners(2,:),'k-x')
    plot(k(1,in,1),k(2,in,1),'rx')
    plot(k(1,on,1),k(2,on,1),'r^')
    plot(k(1,symm,1),k(2,symm,1),'rs')
    set(gca,'FontSize',16)
    title('k-mesh nach Monkhorst-Pack in einem Sechstel der BZ')
    xlabel('k_x')
    ylabel('k_y')
    legend({'Reduzierte 1. BZ',...
        'Gewicht: 6', 'Gewicht: 3', 'Gewicht: 1'},'location','northeast')
    set(gcf, 'Color', 'w');
    axis([min(k(1,:,1)) - 0.1, max(k(1,:,1)) + 0.1, ...
        min(k(2,:,1)) - 0.1, max(k(2,:,1)) + 0.1])
    axis equal
    
    if Ctrl.plot.save == 1
        export_fig(['/home/pmirbach/Masterarbeit/MoS2/Anlauf2/', ...
            'Ergebnisse/k_mesh_red'], '-png');
    end
end

if Ctrl.plot.k_mesh(2) == 1     % Plot über BZ mit Gewichten und Indizierung
    corners = [Para.BZred.symmpts{2} ...
        * [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], [0; 0]];
    [corners] = red_to_BZ(corners);
        
    figure
    hold on
    set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    
    colors = {'b','r','y','g','c','m'};
    marker = {'h','p','x','+','^','v'};
    
    for ii = 1:6
        in = wk == 6;
        on = wk == 3;
        symm = wk == 1;
        
        instr = strcat(colors{ii},marker{1+mod(ii,2)});
        onstr = strcat(colors{ii},marker{3+mod(ii,2)});
        symmstr = strcat(colors{ii},marker{5+mod(ii,2)});
        
        plot(corners(1,:,ii),corners(2,:,ii),'color',0.9 * [1 1 1])
        
        plot(k(1,in,ii),k(2,in,ii),instr)
        plot(k(1,on,ii),k(2,on,ii),onstr)
        plot(k(1,symm,ii),k(2,symm,ii),symmstr)
        
    end
    axis([min(min(k(1,:,:))) - 0.2, max(max(k(1,:,:))) + 0.2, ...
        min(min(k(2,:,:))) - 0.2, max(max(k(2,:,:))) + 0.2])
    axis equal
    set(gca,'FontSize',16)
    title('k-mesh in BZ mit markierten äquivalenten k-Punkten')
    xlabel('k_x')
    ylabel('k_y')
    
    set(gcf, 'Color', 'w');
    
    ind = round(rand(1,4) * size(k,2));
    for jj = 1:numel(ind)
        for ii = 1:6
            plot(k(1,ind(jj),ii),k(2,ind(jj),ii),'ko','markers',8)
            text(k(1,ind(jj),ii) + 0.2, k(2,ind(jj),ii) + 0.2, num2str(jj))
        end
    end
    
    if Ctrl.plot.save == 1
        export_fig(['/home/pmirbach/Masterarbeit/MoS2/Anlauf2/', ...
            'Ergebnisse/k_mesh_index'], '-png');
    end
end