function [Ek, Ev, EGap, Ev_noSOC, H_grad_kx, H_grad_ky, Ek_old] = tight_binding_roesner_1d(Ctrl, Para, k, W90Data)

% D = [cos( Para.k.alpha ) -sin( Para.k.alpha ) ; sin( Para.k.alpha ) cos( Para.k.alpha )];
% for ni = 1:6
%     k(:,:,ni) = D * k(:,:,ni);
% end

TM = W90Data.kLat(1:2,1:2).' * 10;


% Vorläufiger SOC Ansatz:
lambda = Ctrl.material_lambda;                                     % ???  Daniel: 0.074 ??
L_z = -[0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = lambda / 2 * L_z * Para.energy_conversion;


% Berechnung der Eigenwerte und Eigenvektoren ueber das k-mesh
Ek = zeros(6, Para.nr.k);
Ev = zeros(6, 6, Para.nr.k);

% Berechnung der Tight Binding Bandstruktur ohne SOC für die Berechnung der
% Coulomb Matrix Elemente:
Ek_noSOC = zeros(3, Para.nr.k);
Ev_noSOC = zeros(3, 3, Para.nr.k);



HH_TB = zeros(Para.nr.k ,3, 3);

k_m = ( TM \ k(:,:) ).' ;         % Basiswechsel in Vielfache von G
    
[HH_TB(:,:,:), H_grad_kx, H_grad_ky] = getW90Hamiltonian(W90Data, k_m );     %get hamiltonian for every kk point



H_grad_kx = permute( reshape( H_grad_kx,Para.nr.k,9 ),[2,1] ) * Para.energy_conversion;
H_grad_ky = permute( reshape( H_grad_ky,Para.nr.k,9 ),[2,1] ) * Para.energy_conversion;

HH_TB = permute( HH_TB, [2 3 1] );


for nk = 1:Para.nr.k
    
    [ Ek_noSOC(:,nk) , Ev_noSOC(:,:,nk) ] = solve_sort_eig( HH_TB(:,:,nk) );
    
    [ Ek(1:3,nk) , Ev(1:3,1:3,nk) ] = solve_sort_eig( HH_TB(:,:,nk) + H_SOC );
    [ Ek(4:6,nk) , Ev(4:6,4:6,nk) ] = solve_sort_eig( HH_TB(:,:,nk) - H_SOC );
    
end
        
Ek_old = Ek - max( max( Ek( 1, : ) ) );

% Check band gap
% EGap_noSOC = min( Ek_noSOC( 2, : ) - Ek_noSOC( 1, : ) );

EGap = min( min( Ek( [2 5], : ) - Ek( [1 4], : ) ) ) + W90Data.EGapCorr * Para.energy_conversion;

Ek( Para.TB_ind{1}, : ) = - Ek( Para.TB_ind{1}, : ) + max( max( Ek( Para.TB_ind{1}, : ) ) );
Ek( Para.TB_ind{2}, : ) = Ek( Para.TB_ind{2}, : ) - min( min( Ek( Para.TB_ind{2}, 1 ) ) );

