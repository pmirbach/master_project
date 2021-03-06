function handle_path = plot_path(Ctrl,Para,k,Daten,nrpts)

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')

b1 = [Para.k.GV(:,1)];
b2 = [Para.k.GV(:,2)];

k_BZ = reshape(k,2,[]);
k_interp = [k_BZ, k_BZ + b1 * ones(1,size(k_BZ,2)),...
    k_BZ + b2 * ones(1,size(k_BZ,2)),...
    k_BZ + (b1 + b2) * ones(1,size(k_BZ,2))];

[k_p, path_tick] = k_path(Para.BZred.symmpts, Ctrl.plot.path, nrpts);
k_p_tick = 0:size(k_p,2)-1;

handle_path = figure;
hold on
box on
handle_path.Units = 'normalized';
handle_path.Position = [0.1 0.1 0.8 0.8];
handle_path.Color = [1 1 1];

% colors = {'b','r','k','g','c','m'};
colors = {'k','k','k','r','r','r'};
colors = [colors, colors];

linest = {'-','-','-','-','-','-','--','--','--','--','--','--'};

for ii = 1:size(Daten,1)
    Daten_interp = repmat(Daten(ii,:),1,6*4);
    f = scatteredInterpolant(k_interp(1,:)',k_interp(2,:)',...
        Daten_interp(:),'linear','none');
    plot(k_p_tick,f(k_p'),'color',colors{ii},'linestyle',linest{ii},'LineWidth',2)
end
ax = gca;
ax.FontSize = 18;

ax.XTick = path_tick;
ax.XTickLabel = Ctrl.plot.path;

dd = ( max(Daten(:)) - min(Daten(:)) ) / 10;
axy_limits = [min(Daten(:)) - dd, max(Daten(:)) + dd];
ylim(axy_limits)
xlim([min(k_p_tick), max(k_p_tick)]);

for ii = 1:numel(Ctrl.plot.path)
    plot([path_tick(ii), path_tick(ii)],axy_limits,'color', [0.5 0.5 0.5])
end

warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')