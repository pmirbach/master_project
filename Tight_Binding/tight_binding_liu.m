function [Ek, Ev, Ek_noSOC, Ev_noSOC, H_grad_kx , H_grad_ky] = tight_binding_liu(Ctrl, Para, Data)


% Vorläufiger SOC Ansatz:
lambda = Para.TB(end);
L_z = [0 0 0; 0 0 2i; 0 -2i 0];
H_SOC = lambda / 2 * L_z * Para.energy_conversion;


% Berechnung der Eigenwerte und Eigenvektoren über das k-mesh
Ek = zeros(6, Para.nr.k, 6);
Ev = zeros(6, 6, Para.nr.k, 6);

% Berechnung der Tight Binding Bandstruktur ohne SOC für die Berechnung der
% Coulomb Matrix Elemente:
Ek_noSOC = zeros(3, Para.nr.k, 6);
Ev_noSOC = zeros(3, 3, Para.nr.k, 6);

HH_TB = zeros(3, 3, Para.nr.k ,6);
for ni = 1:6
    
    HH_TB(:,:,:,ni) = TB_Liu_TNN_fun(Data.k(:,:,ni), Para.TB);
    HH_TB(:,:,:,ni) = Para.energy_conversion * HH_TB(:,:,:,ni);
    HH_TB(abs( HH_TB ) < 1e-8) = 0;
    
    for nk = 1:Para.nr.k
        
        
        H_TB = HH_TB(:,:,nk,ni);
        
        [Ev_noSOC(:,:,nk,ni), D] = eig( H_TB );       
        Ek_noSOC(:,nk,ni) = diag(D);
        
        
        
        H_TB_SOC_up = (H_TB + H_SOC );               
        H_TB_SOC_down = (H_TB - H_SOC );
        
        [Ev(1:3,1:3,nk,ni) , D_up] = eig(H_TB_SOC_up);
        [Ev(4:6,4:6,nk,ni) , D_down] = eig(H_TB_SOC_down);
        
        Ek(1:3,nk,ni) = real(diag(D_up));
        Ek(4:6,nk,ni) = real(diag(D_down));
        
        
    end
    
end

% % Tests
% k_test_ind = [ Para.symm_indices, round(rand(1,7) * ( Para.nr.k - 1 ) ) + 1 ];   % 10 test kpts including high symmetrie points.
% for ii = 1:size( k_test_ind , 2 )
%     tri = round(rand(1,2) * ( 6 - 1 )) + 1;
%     for ni = 1:size( tri , 2 )
%         [H_TB] = TB_Liu_TNN_fun(Data.k(:,nk,ni), Para.TB);
%         
%     end
% end



Ek = Ek - max(Ek(1,:));
Ek_noSOC = Ek_noSOC - max(Ek_noSOC(1,:));

[H_grad_kx , H_grad_ky] = grad_TB_Liu_TNN_fun(Data.k(1:2,:,1),Para.TB);
