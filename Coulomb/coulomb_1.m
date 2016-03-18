function [Ek_hf,Ek_h,Ek_f] = coulomb_1(constAg,Parameter,Data,CV)

% profile on

plottt = 0;
plot_voranoi = 0;

% k = Data.k;

vorf = constAg.ec^2 / ( 2 * constAg.eps_0 * Parameter.area_real) ...
    / ( 4 * pi )^2;

% Zuerst mal für eine feste Kombination aus Bändern:
% Fock-artig zwischen spin up valenz und leitungsband 1
% Data.Ev(:,1,:) und Data.Ev(:,2,:)

%  V_{k k' k k'}^{l l' l l'}   hier l = 1  l' = 2

Ek_h = Data.Ek;
Ek_f = Data.Ek;
Ek_hf = Data.Ek;

renorm_sign = ones(6);
renorm_sign([1,4],[2,3,5,6]) = -1;
renorm_sign([2,3,5,6],[1,4]) = -1;

coul_diad_h = zeros(6,6,6);
coul_diad_f = zeros(6,6,6);

[V_orbital_h] = fun_coul_orbital_hartree(Parameter.coul_screened);



% tic

% for nk = size(Data.k,2)
for nk = 1:size(Data.k,2)
    
    disp(nk)
    
    % Neue k nach Umklapp Prozess
    k_shift = umklapp1(Parameter, Data.k(1:2,:,:), Data.k(1:2,nk,1));
    
    for nks = 1:size(Data.k,2)
                
        q_v = repmat(Data.k(1:2,nk,1),1,1,6) - k_shift(1:2,nks,:);
        q = squeeze(sqrt( q_v(1,1,:).^2 + q_v(2,1,:).^2 ));
        
%         [V_orbital_f] = fun_coul_orbital_fock(q, ...
%             Parameter.coul_screened, Parameter.coul_kappa );
        [V_orbital_f] = fun_coul_orbital_fock2(q, ...
            Parameter.coul_screened, Parameter.coul_kappa);
                      
        for l1 = 1:6            % Zu renormierendes Band
            
            for l2 = 1:6        % Andere Bänder   

%                 [coul_diad_h] = ...
%                     coul_hartree(Data.Ev(:,:,nk,1), Data.Ev(:,:,nks,:), l1, l2);
%                 
                for ntri = 1:6
                    
                    coul_diad_h(:,:,ntri) = diag(CV(:,:,l1,l1,nk,1)) * ...
                        diag(CV(:,:,l2,l2,nks,ntri))';
                    coul_diad_f(:,:,ntri) = CV(:,:,l1,l1,nk,1) .* ...
                        CV(:,:,l2,l2,nks,ntri).';
                    
                end
                
%                 [coul_diad_f] = ...
%                     coul_fock(Data.Ev(:,:,nk,1), Data.Ev(:,:,nks,:), l1, l2);
                
                               
                V_h = real( sum( sum( sum( coul_diad_h .* V_orbital_h ) ) ) ); 
                V_f = real( sum( sum( sum( coul_diad_f .* V_orbital_f ) ) ) );
                          
                Ek_h(l1,nk) = Ek_h(l1,nk) + renorm_sign(l1,l2) * ...
                    Data.k(3,nks,1) * vorf * V_h * Data.fk(l2,nks);
                Ek_f(l1,nk) = Ek_f(l1,nk) + renorm_sign(l1,l2) * ...
                    Data.k(3,nks,1) * vorf * ( - V_f ) * Data.fk(l2,nks);
                Ek_hf(l1,nk) = Ek_hf(l1,nk) + renorm_sign(l1,l2) * ...
                    Data.k(3,nks,1) * vorf * ( V_h - V_f ) * Data.fk(l2,nks);
                       
            end
            
        end  
                
              
    end

end

% toc
% 
% profile viewer
% profile off

1






% Berechnung der Coulomb WW

if plottt == 1
    
    % Irgend ein festes k'kneu
    % Zufälliger k-Vektor:
    k0_ind = round(rand * size(Data.k,2));
    % k0_ind = 115;
    k0 = Data.k(1:2,k0_ind);
    
    
%     k_shift = umklapp1(Parameter,Data.k,k0);
    k_shift = umklapp1(Parameter, Data.k(1:2,:,:), k0);
    
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
        in = Data.k(3,:,ii) == 6 * Parameter.area_sBZ;
        on = Data.k(3,:,ii) == 3 * Parameter.area_sBZ;
        symm = Data.k(3,:,ii) == 1 * Parameter.area_sBZ;
        
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




