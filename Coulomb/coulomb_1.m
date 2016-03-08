function [] = coulomb_1(constAg,Parameter,Data)

plottt = 0;
plot_voranoi = 0;

k = Data.k;

kappa = 0.1;
vorf = constAg.ec^2 / ( 2 * constAg.eps_0 * Parameter.area_real);

% Zuerst mal für eine feste Kombination aus Bändern:
% Fock-artig zwischen spin up valenz und leitungsband 1
% Data.Ev(:,1,:) und Data.Ev(:,2,:)

%  V_{k k' k k'}^{l l' l l'}   hier l = 1  l' = 2

l1 = 1;
l2 = 2;

V_hartree_12 = zeros(size(k,2));
V_fock_12 = zeros(size(k,2));

tic
for nk = 1:size(k,2)
    % Neue k nach Umklapp Prozess
    kneu = umklapp1(Parameter,k,k(:,nk,1));
    
    parfor nks = 1:size(k,2)
        
        V_hartree_ges = 0;
        V_fock_ges = 0;
        
        for ii = 1:6
            q = norm(k(1:2,nk,1) - kneu(1:2,nks,ii));
            [V_hartree, V_fock] = ...
                fun_Coul_matrix(q,Parameter.coul_screened,kappa);
            
            V_hartree = vorf * V_hartree;
            V_fock = vorf * V_fock;
            
            V_hartree_ges = V_hartree_ges + V_hartree;
            V_fock_ges = V_fock_ges + V_fock;
        end
        
        
        
        V_hartree_12(nk,nks) = V_hartree_ges(1,1);
        V_fock_12(nk,nks) = V_fock_ges(1,1);
        
    end
%     1
end
toc


% Irgend ein festes k'kneu
% Zufälliger k-Vektor:
k0_ind = round(rand * size(k,2));
% k0_ind = 115;
k0 = k(:,k0_ind);



% Berechnung der Coulomb WW


if plottt == 1
    kneu = umklapp1(Parameter,k,k0);
    
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




