function [Ek, Ev, Ek_noSOC, Ev_noSOC] = tight_binding_roesner(Ctrl, Para, Data)


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


tmp = complex( zeros( 3 ) );
for ni = 1:6
    tic
    k_m = [ round( b_m \ Data.k(:,:,ni) ).' , zeros(Para.nr.k,1) ];         % Basiswechsel in Vielfache von G / qr
    
    %get hamiltonian for every kk point
    HH = getW90Hamiltonian(W90Data, k_m / Ctrl.k_mesh_mp.qr);
    toc
    
    tic
    for nk = 1:Para.nr.k
        
        tmp(:,:) = HH(nk, :, :);
       
        H_TB_SOC_up = (tmp + H_SOC );         
        H_TB_SOC_down = (tmp - H_SOC );
            
%         [Ev(1:3,1:3,nk,ni) , D_up] = eig(H_TB_SOC_up);
%         [Ev(4:6,4:6,nk,ni) , D_down] = eig(H_TB_SOC_down);
        
        [ef_up, ev_up] = eig( H_TB_SOC_up );
        [ef_down, ev_down] = eig( H_TB_SOC_down );
        
        [ev_up, I] = sort(diag(real(ev_up)));
        ef_up = ef_up(:, I);
        
        [ev_down, I] = sort(diag(real(ev_down)));
        ef_down = ef_down(:, I);
        
        
        Ev(1:3,1:3,nk,ni) = ef_up;
        Ev(4:6,4:6,nk,ni) = ef_down;
        

        Ek(1:3,nk,ni) = ev_up;
        Ek(4:6,nk,ni) = ev_down;
        
%         Ek(1:3,nk,ni) = real(diag(D_up));
%         Ek(4:6,nk,ni) = real(diag(D_down));
        
        
%         [Ev_noSOC(:,:,nk,ni), D] = eig( tmp );
%         Ek_noSOC(:,nk,ni) = diag(D);
        
        [ef, ev] = eig( tmp );
        
        [ev, I] = sort(diag(real(ev)));
        ef = ef(:, I);
        
        Ev_noSOC(:,:,nk,ni) = ef;
        Ek_noSOC(:,nk,ni) = ev;
        
    end
    toc
    
end

Ek = Ek - max(Ek(1,:));