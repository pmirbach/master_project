function [Ek, Ev, EGap, Ev_noSOC, H_grad_kx, H_grad_ky, Ek_old] = tight_binding_roesner(Ctrl, Para, k, W90Data)

D = [cos( Para.k.alpha ) -sin( Para.k.alpha ) ; sin( Para.k.alpha ) cos( Para.k.alpha )];
for ni = 1:6
    k(:,:,ni) = D * k(:,:,ni);
end

TM = W90Data.kLat(1:2,1:2).' * 10;


if Ctrl.TB_SOC_k
    K_symm_rot = D * Para.BZred.symmpts{2}(:,[1,3]);            % K, K' in Malte's BZ
    lambda = get_SOC_lambda( Para.SOC.lambda_0 , k(:,:,1) , K_symm_rot );
else
    lambda = Para.SOC.lambda_0 * ones(1,Para.nr.k);
end
L_z = -[0 0 0; 0 0 1i; 0 -1i 0];


% Berechnung der Eigenwerte und Eigenvektoren ueber das k-mesh
Ek = zeros(6, Para.nr.k, 6);
Ev = zeros(6, 6, Para.nr.k, 6);

% Berechnung der Tight Binding Bandstruktur ohne SOC für die Berechnung der
% Coulomb Matrix Elemente:
Ek_noSOC = zeros(3, Para.nr.k, 6);
Ev_noSOC = zeros(3, 3, Para.nr.k, 6);



if Ctrl.TB_t_symm == 1   
    if ~strcmp( Ctrl.k_mesh.type , 'symm' )
    	error(' Time symmetrization only possible with k-mesh "symm" ')
    end
    tri_max = 3;
else
    tri_max = 6;
end



HH_TB = zeros(Para.nr.k ,3, 3, 6);
for ni = 1:tri_max
    
    k_m = ( TM \ k(:,:,ni) ).' ;         % Basiswechsel in Vielfache von G
    
    if ni == 1 
        [HH_TB(:,:,:,ni), H_grad_kx, H_grad_ky] = getW90Hamiltonian(W90Data, k_m );     %get hamiltonian for every kk point
    else
        HH_TB(:,:,:,ni) = getW90Hamiltonian(W90Data, k_m );
    end
    HH_TB(:,:,:,ni) = HH_TB(:,:,:,ni) * Para.energy_conversion;
end

% Hamiltonian calculated in eV. 
% Lenght in angs. -> [d/dk] = angs. = 1/10 nm
H_grad_kx = permute( reshape( H_grad_kx , Para.nr.k , 9 ) , [2,1] ) * Para.energy_conversion / 10;
H_grad_ky = permute( reshape( H_grad_ky , Para.nr.k , 9 ) , [2,1] ) * Para.energy_conversion / 10;


HH_TB = permute( HH_TB, [2 3 1 4] );

for ni = 1:tri_max

    for nk = 1:Para.nr.k
        
        [ Ek_noSOC(:,nk,ni) , Ev_noSOC(:,:,nk,ni) ] = solve_sort_eig( HH_TB(:,:,nk,ni) );
        
        H_SOC = lambda(nk) * L_z;   % SOC Hamiltonian with lambda(k)
        
        [ Ek(1:3,nk,ni) , Ev(1:3,1:3,nk,ni) ] = solve_sort_eig( HH_TB(:,:,nk,ni) + H_SOC );
        [ Ek(4:6,nk,ni) , Ev(4:6,4:6,nk,ni) ] = solve_sort_eig( HH_TB(:,:,nk,ni) - H_SOC );
        
    end
    
    if Ctrl.TB_t_symm
        
        Ek_noSOC(:,Para.k_ind.dwn,ni) = Ek_noSOC(:,Para.k_ind.up,ni);   % Unnecessary
        Ek_noSOC(:,:,ni+3) = Ek_noSOC(:,:,ni);
        
        Ek(4:6,Para.k_ind.dwn,ni) = Ek(1:3,Para.k_ind.up,ni);           % Equivalence of spin bands in tri
        Ek(1:3,Para.k_ind.dwn,ni) = Ek(4:6,Para.k_ind.up,ni);    
        Ek(:,Para.k_ind.mid,ni) = ...                                                           % Mean value at middle of redBZ
            repmat( ( Ek(1:3,Para.k_ind.mid,ni) + Ek(4:6,Para.k_ind.mid,ni) ) / 2 , 2 , 1 );    % Unnecessary with H_SOC(k)
        Ek(:,:,ni+3) = Ek(:,:,ni);
                                                                            
        Ev_noSOC(:,:,Para.k_ind.dwn,ni+3) = conj( Ev_noSOC(:,:,Para.k_ind.up,ni) );
        Ev_noSOC(:,:,Para.k_ind.up,ni+3) = conj( Ev_noSOC(:,:,Para.k_ind.dwn,ni) );
        Ev_noSOC(:,:,Para.k_ind.mid,ni+3) = conj( Ev_noSOC(:,:,Para.k_ind.mid,ni) );
        
    end
    
end

%%
% ???? Need rework
% Tests
% sample = get_sample( [ Para.nr.k, Para.nr.tri ] , 10 );
% for ii = 1:size( sample , 1 )
%     nk = sample( ii , 1 );
%     ni = sample( ii , 2 );
%     error_tol = 10 * eps * norm( squeeze( HH_TB(nk,:,:,ni) ) , 2 );
%     ev_test = squeeze( HH_TB(nk,:,:,ni) ) * Ev_noSOC(:,:,nk,ni) - Ev_noSOC(:,:,nk,ni) * diag( Ek_noSOC(:,nk,ni) );
%     if any( ev_test(ev_test>error_tol) )
%         warning('Eigenvectors and Eigenvalues do not solve eigenvalue problem!')
%     end
% end
%%
Ek_old = Ek(:,:,1) - max( max( Ek( Para.TB_ind{1}, : , 1 ) ) );
Ek_old( [2 3 5 6] , : ) = Ek_old( [2 3 5 6] , : ) + W90Data.EGapCorr * Para.energy_conversion;

% Check band gap
EGap_noSOC = min( Ek_noSOC( 2, : , 1 ) - Ek_noSOC( 1, : , 1 ) );
if round( EGap_noSOC , 1 ) ~= round( W90Data.EGapAct * Para.energy_conversion , 1 )
    warning( 'Calculated band gap (%.1f meV) does not agree with Malte (%.1f meV)!' , ...
        round( EGap_noSOC , 1), round( W90Data.EGapAct * Para.energy_conversion , 1 ) )
end

EGap = min( min( Ek( [2 5], : , 1 ) - Ek( [1 4], : , 1 ) ) ) + W90Data.EGapCorr * Para.energy_conversion;

Ek( Para.TB_ind{1}, : , : ) = - Ek( Para.TB_ind{1}, : , : ) + max( max( Ek( Para.TB_ind{1}, : , 1 ) ) );
Ek( Para.TB_ind{2}, : , : ) = Ek( Para.TB_ind{2}, : , : ) - min( min( Ek( Para.TB_ind{2}, Para.k_ind.symm , 1 ) ) );

