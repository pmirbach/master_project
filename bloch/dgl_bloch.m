function dpsik_E = dgl_bloch(t,psik_E,Bloch)

% Übergange: 1 -> 2, 1 -> 3, 4 -> 5, 4 -> 6 

% Berechnung des E-Feldes:
E_t = Bloch.E0 * exp( -1/2 * ( ( t - Bloch.t_peak ) / Bloch.sigma )^2 * 4 * log(2) );


% Extraktion der Psis:
psi_k = psik_E(1:Bloch.nrd * Bloch.nrk);

% Berechnung von        hbar * Omega:
if Bloch.coul_ctrl == 1
    hbarOmega = E_t * conj( Bloch.dipol ) + 1 / ( 2 * pi )^2 * Bloch.coulomb.' * ( psi_k .* Bloch.wkentire ) ;
else 
    hbarOmega = E_t * conj( Bloch.dipol );
end


% Bloch Gleichung:
% dpsik = -1i / Bloch.hbar * Bloch.Esum .* psi_k ...
%     + 1i * hbarOmega / Bloch.hbar - Bloch.gamma / Bloch.hbar * psi_k;

dpsik = -1i / Bloch.hbar * ( Bloch.Esum .* psi_k - ( hbarOmega ) - 1i * Bloch.gamma * psi_k );



% Berechnung der gesamten Polarisation:
P_t = 1 / (2 * pi)^2 * sum( conj( Bloch.dipol(:,1) ) .* psi_k .* Bloch.wk );

% Berechnung der Fourier Transformation:
fw = exp(1i * Bloch.w * t);

d_P_w = P_t * fw;
d_E_w = E_t * fw;

% Zeitentwicklung:
dpsik_E = [dpsik; d_P_w; d_E_w];
