function [Ek, Ev, Ek_noSOC, Ev_noSOC, HH_grad_kx, HH_grad_ky] = tight_binding_roesner(Ctrl, Para, Data)


a0_form = num2str( 10 * Ctrl.lattice_constant , '%.3f' );
seed = [ './Tight_Binding/ab_initio/02_Materials/', Ctrl.material, '/a0_', a0_form, 'A/01_WannierTB/02_G0W0/02_Mo3d/wannier90' ];

SOCSettings.type  = 'none'; % 'first' or 'second' or 'none'

% load wannier90 Data
W90Data = loadW90Data(seed, SOCSettings);


% Vorläufiger SOC Ansatz:
% SOC Hamiltonian
lambda = Para.TB(end);
L_z = [0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = lambda / 2 * L_z;


% Berechnung der Eigenwerte und Eigenvektoren ueber das k-mesh
Ek = zeros(6, Para.nr.k, 6);
Ev = zeros(6, 6, Para.nr.k, 6);

% Berechnung der Tight Binding Bandstruktur ohne SOC für die Berechnung der
% Coulomb Matrix Elemente:
Ek_noSOC = zeros(3, Para.nr.k, 6);
Ev_noSOC = zeros(3, 3, Para.nr.k, 6);


% Basiswechsel in Vielfache von den reziproken Gittervektoren:
b_m = abs( Para.k.GV ) / Ctrl.k_mesh_mp.qr;    % Maltes GV / qr


H_TB = complex( zeros( 3 ) );
for ni = 1:6
    
    k_m = [ round( b_m \ Data.k(:,:,ni) ).' , zeros(Para.nr.k,1) ];         % Basiswechsel in Vielfache von G / qr
    
    if ni == 1
        %get hamiltonian for every kk point
        [HH, HH_grad_kx, HH_grad_ky] = getW90Hamiltonian(W90Data, k_m / Ctrl.k_mesh_mp.qr);
    else
        HH = getW90Hamiltonian(W90Data, k_m / Ctrl.k_mesh_mp.qr);
    end
    
    for nk = 1:Para.nr.k
        
        H_TB(:,:) = HH(nk, :, :);
        
        [ef, ev] = eig( H_TB *1e3 );
        
        [Ek_noSOC(:,nk,ni), I] = sort(diag(real(ev)));
        Ev_noSOC(:,:,nk,ni) = ef(:, I);
        
        
        H_TB_SOC_up = (H_TB + H_SOC ) *1e3;
        H_TB_SOC_down = (H_TB - H_SOC ) *1e3;
        
        
        [ef_up, ev_up] = eig( H_TB_SOC_up );
        [ef_down, ev_down] = eig( H_TB_SOC_down );
        
        [Ek(1:3,nk,ni), I] = sort(diag(real(ev_up)));
        Ev(1:3,1:3,nk,ni) = ef_up(:, I);
        
        [Ek(4:6,nk,ni), I] = sort(diag(real(ev_down)));
        Ev(4:6,4:6,nk,ni) = ef_down(:, I);
        
    end
    
end


% Ev = Ev( [1,3,2,6,5,4],:,:,: );


Ek = Ek - max(Ek(1,:));