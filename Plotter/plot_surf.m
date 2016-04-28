function handle_surf = plot_surf(Para,k,Daten,layout)

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')

% Skalare und Vektoren zur Erzeugung der Subplots
Nplots = size(Daten,1);
Nsubplots = layout(1) * layout(end);
subii = repmat(1:Nsubplots,1,ceil(Nplots / Nsubplots));
Nr_fig = ceil(Nplots / Nsubplots);
x2 = ones(Nsubplots,1)*(1:Nr_fig);
er = reshape(x2,1,Nr_fig * Nsubplots);

% Reihenfolge der Symm Pts zum Plotten
corners = [Para.BZred.symmpts{2} ...
    * [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], [0; 0]];

if size(k,3) > 1
    x_lim = [ min(min(k(1,:,:)))  , max(max(k(1,:,:)))  ];
    y_lim = [ min(min(k(2,:,:)))  , max(max(k(2,:,:)))  ];
    gen = 300;
    k_plot = reshape(k,3,[]);
    
    [corners] = red_to_BZ(corners);
    corners = reshape(corners,2,[]);
else
    x_lim = [ min(k(1,:,1)) , max(k(1,:,1)) ];
    y_lim = [ min(k(2,:,1)) , max(k(2,:,1)) ];
    gen = 200;
    k_plot = k(:,:,1);
end

% Erzeugung: Meshgrid zum Plotten der Interpolation
xlin = linspace(x_lim(1),x_lim(2),gen);
ylin = linspace(y_lim(1),y_lim(2),gen);
[X,Y] = meshgrid(xlin,ylin);

handle_surf = gobjects(1,Nr_fig);

for ii = 1:Nr_fig
    handle_surf(ii) = figure;
    set(handle_surf(ii), ...
        'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    set(handle_surf(ii), 'Color', 'w');
end


for ii = 1:Nplots
    % Auswahl figure
    figure(handle_surf(er(ii)))
    % Auswahl subplot
    subplot(layout(1),layout(end),subii(ii))
    hold on
    
    if size(k,3) > 1
        Daten_plot = repmat(Daten(ii,:),1,6);
    else
        Daten_plot = Daten(ii,:);
    end
        
    f = scatteredInterpolant(k_plot(1,:)',k_plot(2,:)',...
        Daten_plot(:),'linear','none');
    Z = f(X,Y);
    
    % figure
    over_data = max(max(Daten)) + 10;
    plot3(corners(1,:),corners(2,:),over_data * ones(1,size(corners,2)),'k-x')
    surf(X,Y,Z,'EdgeColor','none','LineStyle','none',...
        'FaceLighting','phong');
    view(2)
    colormap jet
    colorbar
    
    d = [(x_lim(2) - x_lim(1)) / 15, (y_lim(2) - y_lim(1)) / 15];
    axis([x_lim(1)-d(1), x_lim(2)+d(1), y_lim(1)-d(2), y_lim(2)+d(2)]);       
    axis equal
    axis off
    
    text(Para.BZred.symmpts{2}(1,:) + 0.8 * [-0.3 -0.3 -0.3 0.5], ...
        Para.BZred.symmpts{2}(2,:) + 0.8 * [-1 -1 1 0.5], ...
        Para.BZred.symmpts{1})
end
warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')