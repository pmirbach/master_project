function [Ek, Ev] = tight_binding_liu(Ctrl, Parameter, Data)

% Berechnung der Eigenwerte und Eigenvektoren über das k-mesh
Ek = zeros(6, Parameter.nrk, 6);
Ev = zeros(6, 6, Parameter.nrk, 6);

% SOC Hamiltonian
lambda = Parameter.TB.liu.values(end);
L_z = [0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = Ctrl.SOC * lambda / 2 * L_z;
% H_SOC = Ctrl.SOC * lambda / 2 * blkdiag(L_z, -L_z);

for nk = 1:Parameter.nrk
    
    for ni = 1:6
        
        switch Ctrl.method
            case 'NN'
                [H_TB] = TB_Liu_NN_fun(Data.k(1:2,nk,ni), Parameter);
            case 'TNN'
                [H_TB] = TB_Liu_TNN_fun(Data.k(1:2,nk,ni), Parameter);
        end

        %     H_TB_SOC = blkdiag(H_TB, H_TB) + H_SOC;
        H_TB_SOC_up = (H_TB + H_SOC) * 1e3;         % Arbeiten in meV
        H_TB_SOC_down = (H_TB - H_SOC) * 1e3;
                
        %     [Ev(:,:,ii) , D] = eig(H_TB_SOC);
        [Ev(1:3,1:3,nk,ni) , D_up] = eig(H_TB_SOC_up);
        [Ev(4:6,4:6,nk,ni) , D_down] = eig(H_TB_SOC_down);
        
        %     Ek(:,ii) = diag(D);
        Ek(1:3,nk,ni) = diag(D_up);
        Ek(4:6,nk,ni) = diag(D_down);
        
    end
    
end

Ek = Ek - max(Ek(1,:));