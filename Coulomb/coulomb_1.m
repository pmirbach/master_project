function [] = coulomb_1(Parameter,k)

plottt = 1;
plot_voranoi = 0;

% size(k);

% Irgend ein festes k'
% Zufälliger k-Vektor:

k0_ind = round(rand * size(k,2));
% k0_ind = 115;
k0 = k(:,k0_ind);

kneu = umklapp1(Parameter,k,k0);


if plottt == 1
    
    corners = Parameter.symmpts{2}(:,[2 3]);
    [corners] = red_to_BZ(corners);
    corners = reshape(corners,2,[]);
    corners_round = round(corners*1e10)*1e-10;
    [dummy1, m2, dummy2] = unique(corners_round','rows');
    corners = [corners(:,sort(m2)), corners(:,1)];
    
    figure
    hold on
    set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    set(gcf,'Color','w')
    plot(corners(1,:),corners(2,:),'color',0.7 * [0 1 1],'linewidth',2);
    axis equal
    plot(corners(1,:),corners(2,:));
    
    colors = {'b','r','y','g','c','m'};
    marker = {'h','p','x','+','^','v'};
    
    for ii = 1:6
        in = kneu(3,:,ii) == 6 * Parameter.area_sBZ;
        on = kneu(3,:,ii) == 3 * Parameter.area_sBZ;
        symm = kneu(3,:,ii) == 1 * Parameter.area_sBZ;
        
        instr = strcat(colors{ii},marker{1+mod(ii,2)});
        onstr = strcat(colors{ii},marker{3+mod(ii,2)});
        symmstr = strcat(colors{ii},marker{5+mod(ii,2)});
        
        %     plot(corners(1,:,ii),corners(2,:,ii),'color',0.9 * [1 1 1])
        
        plot(kneu(1,in,ii),kneu(2,in,ii),instr)
        plot(kneu(1,on,ii),kneu(2,on,ii),onstr)
        plot(kneu(1,symm,ii),kneu(2,symm,ii),symmstr)
        
    end
    plot(k0(1),k0(2),'ko')
    corners0 = corners + k0(1:2)*ones(1,size(corners,2));
    plot(corners0(1,:),corners0(2,:),'color',0.4 * [1 1 1],'linewidth',1)
end

if plot_voranoi == 1
    figure;
    voronoi(k(1,:),k(2,:));
    hold on
    axis([min(k(1,:))-0.1, max(k(1,:))+0.1, min(k(2,:))-0.1, max(k(2,:))+0.1])
    corners = [Parameter.symmpts{2} ...
        * [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], [0; 0]];
    plot(corners(1,:),corners(2,:),'k-x')
end