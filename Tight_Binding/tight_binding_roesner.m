function [Ek, Ev, Ek_noSOC, Ev_noSOC, H_grad_kx, H_grad_ky] = tight_binding_roesner(Ctrl, Para, Data, W90Data)


% Vorläufiger SOC Ansatz:
lambda = 0.073;                                     % ???  Daniel: 0.074 ??
L_z = [0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = lambda / 2 * L_z * Para.energy_conversion;


% Berechnung der Eigenwerte und Eigenvektoren ueber das k-mesh
Ek = zeros(6, Para.nr.k, 6);
Ev = zeros(6, 6, Para.nr.k, 6);

% Berechnung der Tight Binding Bandstruktur ohne SOC für die Berechnung der
% Coulomb Matrix Elemente:
Ek_noSOC = zeros(3, Para.nr.k, 6);
Ev_noSOC = zeros(3, 3, Para.nr.k, 6);


k_m = zeros( Para.nr.k , 3 , 6 );

if ~Ctrl.cmp.use_k
    % Basiswechsel in Vielfache von den reziproken Gittervektoren:
    b_m = abs( Para.k.GV ) / Ctrl.k_mesh_mp.qr;                                 % Maltes GV / qr
    for ni = 1:6
        k_m(:,1:2,ni) = round( b_m \ Data.k(:,:,ni) ).' / Ctrl.k_mesh_mp.qr;    % Basiswechsel in Vielfache von G (hohe Genauigkeit)
    end
elseif Ctrl.cmp.use_k
    b_m = abs( Para.k.GV );                         % Maltes GV 
    for ni = 1:6
        k_m(:,1:2,ni) = ( b_m \ Data.k(:,:,ni) ).' ;         % Basiswechsel in Vielfache von G
    end
else
    error('Something strange with kpts compare in TB_ab_initio')
end


HH_TB = zeros(Para.nr.k ,3, 3, 6);
H_TB = complex( zeros( 3 ) );
for ni = 1:6
        
    if ni == 1
        %get hamiltonian for every kk point
        [HH_TB(:,:,:,ni), H_grad_kx, H_grad_ky] = getW90Hamiltonian(W90Data, k_m(:,:,ni) );
    else
        HH_TB(:,:,:,ni) = getW90Hamiltonian(W90Data, k_m(:,:,ni) );
    end
    HH_TB(:,:,:,ni) = HH_TB(:,:,:,ni) * Para.energy_conversion;
    
    for nk = 1:Para.nr.k
        
        H_TB(:,:) = HH_TB(nk,:,:,ni);
        
        [ef, ev] = eig( H_TB );
        
        [Ek_noSOC(:,nk,ni), I] = sort(diag(real(ev)));
        Ev_noSOC(:,:,nk,ni) = ef(:, I);
        
        
        H_TB_SOC_up = (H_TB + H_SOC );
        H_TB_SOC_down = (H_TB - H_SOC );
        
        
        [ef_up, ev_up] = eig( H_TB_SOC_up );
        [ef_down, ev_down] = eig( H_TB_SOC_down );
        
        [Ek(1:3,nk,ni), I] = sort(diag(real(ev_up)));
        Ev(1:3,1:3,nk,ni) = ef_up(:, I);
        
        [Ek(4:6,nk,ni), I] = sort(diag(real(ev_down)));
        Ev(4:6,4:6,nk,ni) = ef_down(:, I);
        
    end
    
end

% Tests
sample = get_sample( [ Para.nr.k, Para.nr.tri ] , 10 );
for ii = 1:size( sample , 1 )
    nk = sample( ii , 1 );
    ni = sample( ii , 2 );
    error_tol = 10 * eps * norm( squeeze( HH_TB(nk,:,:,ni) ) , 2 );
    ev_test = squeeze( HH_TB(nk,:,:,ni) ) * Ev_noSOC(:,:,nk,ni) - Ev_noSOC(:,:,nk,ni) * diag( Ek_noSOC(:,nk,ni) );
    if any( ev_test(ev_test>error_tol) )
        warning('Eigenvectors or Eigenvalues poorly calculated!')
    end
end



% Ev = Ev( [1,3,2,6,5,4],:,:,: );

Ek = Ek - max(Ek(1,:));
Ek_noSOC = Ek_noSOC - max(Ek_noSOC(1,:));





H_grad_kx = permute(reshape(H_grad_kx,Para.nr.k,9),[2,1]) * 1e3;
H_grad_ky = permute(reshape(H_grad_ky,Para.nr.k,9),[2,1]) * 1e3;
