function V_ab_interpl = get_background_screening( Ctrl , Para , Coul_ME , minq )

% Define q for coulomb-interaction
% Find q (unique in minq)
% q_pts = sort( unique( round( minq,12 ) ) );
% q_pts( q_pts == 0 ) = [];
% q = q_pts.';

%% Alternative declaration as a linear q mesh
q_max = max(minq(:));
q_min = min(minq(:));
del = 5000;
q = q_min : (q_max - q_min) / del : q_max;
q( q == 0 ) = [];
%%
% q = linspace(0,27,1000);
% q(1) = [];

%%
nr_q = numel( q );
A = Para.real.area;

% % Memory allocation for U and eps ind diag basis
U_diag = zeros(3,3,nr_q);
eps_diag = zeros(3,3,nr_q);
% % Memory allocation fur U, eps, V in ab base
U_ab_q = zeros(3,3,nr_q);
V_ab_q = zeros(3,3,nr_q);


% % Calculation of U in diag basis
% Quadratic fit for leading U
% U_diag(1,1,:) = 3 * 1 ./ q .* 1 ./ ( 1 + Coul_ME.U.gamma_quad .* q );
% Cubic fit for leading U
U_diag(1,1,:) = 3 * 1 ./ q .* 1 ./ ( 1 + Coul_ME.U.gamma .* q + Coul_ME.U.delta .* q.^2 );

U_diag(2,2,:) = repmat( Coul_ME.U.micro(1) / Para.vorf.coul * A , 1 , nr_q );
U_diag(3,3,:) = repmat( Coul_ME.U.micro(2) / Para.vorf.coul * A , 1 , nr_q );

% % Calculation of eps in diag basis
% Calculate eps_inf and height either with Resta fit or simple 
if Ctrl.Coul.Resta_fit_eps                                                                                                     % Inserted NOT
    
    a = Coul_ME.eps.Resta(1);
    b = Coul_ME.eps.Resta(2);
    c = Coul_ME.eps.Resta(3);
    d = Coul_ME.eps.Resta(4);
    e = Coul_ME.eps.Resta(5);
    
    eps_inf_q = ( a + q.^2 ) ./ ( a / b  * sin( q * c ) ./ ( q * c ) + q.^2 ) + e;
    height = d;
       
else
    
    eps_inf_q = repmat( Coul_ME.eps.inf , 1 , nr_q );
    height = Coul_ME.eps.L_2d_macro;
    
end

% % Old formula without dielectrica
% eps_diag(1,1,:) = eps_inf_q .* ( eps_inf_q + 1 - ( eps_inf_q - 1) .* exp( -q * height ) ) ...
%     ./ ( eps_inf_q + 1 + ( eps_inf_q - 1) .* exp( -q * height ) );

% New formula for leding eps. Includes different dielectrica above and below
g_2 = ( eps_inf_q - Ctrl.Coul.eps_2 ) ./ ( eps_inf_q + Ctrl.Coul.eps_2 );
g_3 = ( eps_inf_q - Ctrl.Coul.eps_3 ) ./ ( eps_inf_q + Ctrl.Coul.eps_3 );
eps_diag(1,1,:) = eps_inf_q .* ( 1 - g_2 .* g_3 .* exp( -2 * q * height ) ) ...
    ./ ( 1 + ( g_2 + g_3 ) .* exp( -q * height ) + g_2 .* g_3 .* exp( -2 * q * height ) );

eps_diag(2,2,:) = repmat( Coul_ME.eps.micro(1) , 1 , nr_q );
eps_diag(3,3,:) = repmat( Coul_ME.eps.micro(2) , 1 , nr_q );



% Ev = Coul_ME.U_Ev ;

Ev2 = Coul_ME.U_Ev;
% Ev2(:,1) = 1/sqrt(3);
% Ev2(:,3) = [ 0; -1/sqrt(2); 1/sqrt(2) ];

Ev = Ev2;

% Ev = Coul_ME.U_Ev([1,3,2],:);

for ii = 1:nr_q
    
    V_diag = U_diag(:,:,ii) / eps_diag(:,:,ii);

%     U_ab_q(:,:,ii) = Ev * U_diag(:,:,ii) * transpose( Ev );
    U_ab_q(:,:,ii) = Ev * ( U_diag(:,:,ii) / Ev );
    
    V_ab_q(:,:,ii) = Ev * V_diag * transpose( Ev );
%     V_ab_q(:,:,ii) = Ev * ( V_diag / Ev );
    
end

U_ab_q = reshape( U_ab_q , 9 , []);
V_ab_q = reshape( V_ab_q , 9 , []);


% V_ab_q( V_ab_q < 0 ) = 0;

V_ab_q( [4 7 8], : ) = [];
U_ab_q( [4 7 8], : ) = [];

% %%%%%%%%%%%%%%%%%%
% Danielq = load('qgrid_631.mat');
% DanielV = load('V_ab_q_von631aus.mat');
% 
% q_d = Danielq.q(631,:,1);
% V_ab_qq = DanielV.V_ab_q;
% 
% V_ab_qq = permute(V_ab_qq,[2,3,1]);
% V_ab_q_d = reshape(V_ab_qq, 9 , []);
% V_ab_q_d( [4 7 8], : ) = [];
% for ii = 1:6
%     subplot(2,3,ii)
% %     plot(q,Para.vorf.coul * V_ab_q(ii,:) - V_ab_q_d(ii,2:end))
%     
%     plot(q,Para.vorf.coul * V_ab_q(ii,:),'-x')
%     hold on
%     plot(q_d,V_ab_q_d(ii,:)./q_d,'x')
% end
% %%%%%%%%%%%%%%%%%%

q = [0 , q];
V_ab_q = [Para.coul.pol * ones(6,1) , V_ab_q];

V_ab_interpl = cell(1,6);
for ii = 1:6
    V_ab_interpl{ii} = griddedInterpolant( q , V_ab_q(ii,:) );    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% F�r die Arbeit!

% titlestring = {'d_{z^2} d_{z^2}','d_{z^2} d_{x^2-y^2}','d_{z^2} d_{xy}','d_{x^2-y^2} d_{x^2-y^2}','d_{x^2-y^2} d_{xy}','d_{xy} d_{xy}'};
% 
% figure(6)
% set(gcf,'color','w');

% for ii = 1:6
%     subplot(2,3,ii)
%     title(titlestring{ii})
%     
%     hold on
%     plot(q,Para.vorf.coul*V_ab_q(ii,:)/Para.energy_conversion)
%     
%     axis([-inf 5 0 20])  
%     
%     ylabel('V in ev')
%     
% %     legend('MoS_2','MoSe_2','WS_2','WSe_2')
% end

% for ii = 1:6
%     subplot(2,3,ii)
%     title(titlestring{ii})
%     
%     hold on
%     plot(q(2:end),U_ab_q(ii,:) ./ V_ab_q(ii,2:end))
%     
%     ylabel(char(949))
% end

% figure(6)
% set(gcf,'units','normalized','position',[.1 .1 .4 .4])
% hold on
% set(gcf,'color','w');
% plot(q(2:end),U_ab_q(1,:) ./ V_ab_q(1,2:end),'k','linewidth',2)
% axis([0 13.17 0 9.5])
% set(gca,'fontsize',18);
% ylabel(char(949))
% xlabel('q in nm^{-1}')
% 
% figure(7)
% set(gcf,'color','w');
% plot(q,Para.vorf.coul*V_ab_q(ii,:)/Para.energy_conversion,'k','linewidth',2)
% set(gca,'fontsize',18);
% axis([-inf 10 0 20])
% ylabel('V in eV')
% xlabel('q in nm^{-1}')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % Tests:

% % Old model data
% V_ab_q_old_model = zeros(6,nr_q+1);
% U_ab_q_old_model = zeros(6,nr_q+1);
% eps_ab_q_old_model = zeros(6,nr_q+1);
% for ii = 1:6
%     [ V_tmp , U_tmp , eps_tmp ] = final_coul_scr( q , Para.coul.screened(ii,:) , Para.coul.pol );
%     V_ab_q_old_model(ii,:) = V_tmp.';
%     U_ab_q_old_model(ii,:) = U_tmp.';
%     eps_ab_q_old_model(ii,:) = eps_tmp.';
% end


% % Compare U
% figure; plot(q(2:end),squeeze(U_diag(1,1,:)))
% title('U_{diag}(1,1)(q)')
% 
% figure
% for ii = 1:6
%     subplot(2,3,ii)
%     plot(q(2:end), U_ab_q(ii,:),'-' )
%     hold on
%     plot(q, U_ab_q_old_model(ii,:), 'r--')
% end
% mtit( 'U_{a,b}(q)' , 'fontsize' , 18 )

% % Compare eps
% figure; plot(q(2:end),squeeze(eps_diag(1,1,:)))
% title('eps_{diag}(1,1)(q)')
% 
% figure
% for ii = 1:6
%     subplot(2,3,ii)
%     plot(q(2:end), eps_ab_q(ii,:),'-' )
%     hold on
%     plot(q, eps_ab_q_old_model(ii,:), 'r--')
% end
% mtit( 'eps_{a,b}(q)' , 'fontsize' , 18 )

% % Compare V
% figure
% for ii = 1:6
%     subplot(2,3,ii)
%     plot(q, V_ab_q(ii,:),'-' )
%     hold on
%     plot(q, V_ab_q_old_model(ii,:), 'r--')
% end


% figure
% for ii = 1:6
%     subplot(2,3,ii)
%     plot(q, V_ab_interpl{ii}(q) )
% end
