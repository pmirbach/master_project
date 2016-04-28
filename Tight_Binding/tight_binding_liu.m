function [Ek, Ev, Ek_noSOC, Ev_noSOC] = tight_binding_liu(Ctrl, Para, Data)

% Berechnung der Eigenwerte und Eigenvektoren über das k-mesh
Ek = zeros(6, Para.nr.k, 6);
Ev = zeros(6, 6, Para.nr.k, 6);

% SOC Hamiltonian
lambda = Para.TB(end);
L_z = [0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = Ctrl.SOC * lambda / 2 * L_z;

% Berechnung der Tight Binding Bandstruktur ohne SOC für die Berechnung der
% Coulomb Matrix Elemente:

Ek_noSOC = zeros(3, Para.nr.k, 6);
Ev_noSOC = zeros(3, 3, Para.nr.k, 6);

for nk = 1:Para.nr.k
    
    for ni = 1:6
        
        [H_TB] = TB_Liu_TNN_fun(Data.k(1:2,nk,ni), Para.TB);
        
        H_TB_SOC_up = (H_TB + H_SOC) * 1e3;         % Arbeiten in meV
        H_TB_SOC_down = (H_TB - H_SOC) * 1e3;
        
        [Ev(1:3,1:3,nk,ni) , D_up] = eig(H_TB_SOC_up);
        [Ev(4:6,4:6,nk,ni) , D_down] = eig(H_TB_SOC_down);
        
        Ek(1:3,nk,ni) = diag(D_up);
        Ek(4:6,nk,ni) = diag(D_down);
        
        [Ev_noSOC(:,:,nk,ni), D] = eig( H_TB * 1e3 );
        Ek_noSOC(:,nk,ni) = diag(D);
        
    end
    
end

Ek = Ek - max(Ek(1,:));