function [Ek, Ev, Ek_noSOC, Ev_noSOC, H_grad_kx, H_grad_ky] = tight_binding_roesner(Ctrl, Para, k, W90Data)

D = [cos( Para.k.alpha ) -sin( Para.k.alpha ) ; sin( Para.k.alpha ) cos( Para.k.alpha )];
for ni = 1:6
    k(:,:,ni) = D * k(:,:,ni);
end

TM = W90Data.kLat(1:2,1:2).' * 10;


% Vorläufiger SOC Ansatz:
lambda = 0.074;                                     % ???  Daniel: 0.074 ??
L_z = -[0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = lambda / 2 * L_z * Para.energy_conversion;


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

HH_TB = permute( HH_TB, [2 3 1 4] );



tri_max = 6;

for ni = 1:tri_max

    for nk = 1:Para.nr.k
        
        [ Ek_noSOC(:,nk,ni) , Ev_noSOC(:,:,nk,ni) ] = solve_sort_eig( HH_TB(:,:,nk,ni) );
        
        [ Ek(1:3,nk,ni) , Ev(1:3,1:3,nk,ni) ] = solve_sort_eig( HH_TB(:,:,nk,ni) + H_SOC );
        [ Ek(4:6,nk,ni) , Ev(4:6,4:6,nk,ni) ] = solve_sort_eig( HH_TB(:,:,nk,ni) - H_SOC );
        
    end
    
    if Ctrl.TB_t_symm == 1
        
%         Ek_noSOC(:,Para.k_ind.dwn,ni) = Ek_noSOC(:,Para.k_ind.up,ni);
%         Ek(4:6,Para.k_ind.dwn,ni) = Ek(1:3,Para.k_ind.up,ni);
%         Ek(1:3,Para.k_ind.dwn,ni) = Ek(4:6,Para.k_ind.up,ni);
        
%         Ek(:,Para.k_ind.mid,ni) = repmat( ( Ek(1:3,Para.k_ind.mid,ni) + Ek(4:6,Para.k_ind.mid,ni) ) / 2 , 2 )
        
    end
    
end



% for ni = 1:tri_max
%             
%     for nk = 1:Para.nr.k
%         
%         H_TB(:,:) = HH_TB(:,:,nk,ni);
%         H_TB_SOC_up = (H_TB + H_SOC );
%         H_TB_SOC_down = (H_TB - H_SOC );
%         
%         
%         [ef, ev] = eig( H_TB );
%         
%         [Ek_noSOC(:,nk,ni), I] = sort(diag(real(ev)));
%         Ev_noSOC(:,:,nk,ni) = ef(:, I);
%         
%         
%         [ef_up, ev_up] = eig( H_TB_SOC_up );
%         [ef_down, ev_down] = eig( H_TB_SOC_down );
%         
%         [Ek(1:3,nk,ni), I] = sort(diag(real(ev_up)));
%         Ev(1:3,1:3,nk,ni) = ef_up(:, I);
%         
%         [Ek(4:6,nk,ni), I] = sort(diag(real(ev_down)));
%         Ev(4:6,4:6,nk,ni) = ef_down(:, I);
%         
%     end
%     
% end





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



Ek = Ek - max(Ek(1,:));
Ek_noSOC = Ek_noSOC - max(Ek_noSOC(1,:));


H_grad_kx = permute(reshape(H_grad_kx,Para.nr.k,9),[2,1]) * 1e3;
H_grad_ky = permute(reshape(H_grad_ky,Para.nr.k,9),[2,1]) * 1e3;
