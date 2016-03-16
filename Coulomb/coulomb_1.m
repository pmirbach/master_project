function [] = coulomb_1(constAg,Parameter,Data)

plottt = 0;
plot_voranoi = 0;

% k = Data.k;

vorf = constAg.ec^2 / ( 2 * constAg.eps_0 * Parameter.area_real);

% Zuerst mal für eine feste Kombination aus Bändern:
% Fock-artig zwischen spin up valenz und leitungsband 1
% Data.Ev(:,1,:) und Data.Ev(:,2,:)

%  V_{k k' k k'}^{l l' l l'}   hier l = 1  l' = 2

Ek_h = Data.Ek;
Ek_f = Data.Ek;
Ek_hf = Data.Ek;


[V_orbital_h] = fun_coul_orbital_hartree(Parameter.coul_screened);

% % % repmat(blkdiag(V_orbital_h,V_orbital_h),1,1,6)

% tic
for nk = 1:size(Data.k,2)
    
%     Data.Ek(:,nk) % ist zu renormieren
    
    % Neue k nach Umklapp Prozess
    k_shift = umklapp1(Parameter, Data.k, Data.k(:,nk,1));
    
    for nks = 2:size(Data.k,2)
                
        q_v = squeeze( repmat(Data.k(1:2,nk,1),1,1,6) - k_shift(1:2,nks,:) );
        q = sqrt( q_v(1,:).^2 + q_v(2,:).^2 );
        
        [V_orbital_f] = fun_coul_orbital_fock(q, ...
            Parameter.coul_screened, Parameter.coul_kappa );
                      
        for l1 = 1:6            % Zu renormierendes Band
            
            for l2 = 1:6        % Andere B�nder
                
%                 [coul_diad_h, coul_diad_f] = ...
%                     test_diad(Data.Ev(:,:,nk), Data.Ev(:,:,nks),[l1, l2, l2, l1]);
                
                [coul_diad_h, coul_diad_f] = ...
                    fun_coul_diad(Data.Ev(:,:,nk,1), Data.Ev(:,:,nks,:), ...
                    [l1, l2, l2, l1]);
                
                V_h = sum(sum(sum(coul_diad_h .* V_orbital_h)))
                V_f = sum(sum(sum(coul_diad_f .* V_orbital_f)))
                
%                 V_h = sum( sum( V_orbital_h .* coul_diad_h ) )
%                 V_f = sum( sum( V_orbital_f .* coul_diad_f ) )
%                 disp(coul_diad_h)
%                 disp(coul_diad_f)

%                 disp(imag(V_h))
%                 disp(imag(coul_diad_f))
%                 disp(V_f)
                
                1;
                
            end
            
        end  
                
              
    end

end
% toc


% Irgend ein festes k'kneu
% Zufälliger k-Vektor:
k0_ind = round(rand * size(Data.k,2));
% k0_ind = 115;
k0 = Data.k(:,k0_ind);



% Berechnung der Coulomb WW

if plottt == 1
    k_shift = umklapp1(Parameter,Data.k,k0);
    
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
        in = k_shift(3,:,ii) == 6 * Parameter.area_sBZ;
        on = k_shift(3,:,ii) == 3 * Parameter.area_sBZ;
        symm = k_shift(3,:,ii) == 1 * Parameter.area_sBZ;
        
        instr = strcat(colors{ii},marker{1+mod(ii,2)});
        onstr = strcat(colors{ii},marker{3+mod(ii,2)});
        symmstr = strcat(colors{ii},marker{5+mod(ii,2)});
        
        %     plot(corners(1,:,ii),corners(2,:,ii),'color',0.9 * [1 1 1])
        
        plot(k_shift(1,in,ii),k_shift(2,in,ii),instr)
        plot(k_shift(1,on,ii),k_shift(2,on,ii),onstr)
        plot(k_shift(1,symm,ii),k_shift(2,symm,ii),symmstr)
        
    end
    plot(k0(1),k0(2),'ko')
    corners0 = corners + k0(1:2)*ones(1,size(corners,2));
    plot(corners0(1,:),corners0(2,:),'color',0.4 * [1 1 1],'linewidth',1)
end

if plot_voranoi == 1
    figure;
    voronoi(Data.k(1,:),Data.k(2,:));
    hold on
    axis([min(Data.k(1,:))-0.1, max(Data.k(1,:))+0.1,...
        min(Data.k(2,:))-0.1, max(Data.k(2,:))+0.1])
    corners = [Parameter.symmpts{2} ...
        * [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], [0; 0]];
    plot(corners(1,:),corners(2,:),'k-x')
end




