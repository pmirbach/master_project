function dpsik_E = dgl_bloch2(t,psik_E,Bloch)


% Extraktion der Psis:
psi_k = psik_E(1:Bloch.nrk);

% Berechnung des E-Feldes:
E_t = Bloch.E0 * exp(-1/2 * ( ( t - Bloch.t_peak ) / Bloch.sigma )^2 * 4 * log(2) );

% Berechnung von        hbar * Omega:
if Bloch.coul_ctrl == 1
    hbarOmega = E_t * conj(Bloch.dipol) + ( 1 / ( 2 * pi )^2 * ( psi_k.' .* Bloch.wkentire ) * Bloch.coulomb.' ).';
else 
    hbarOmega = E_t * conj(Bloch.dipol);
end


% Bloch Gleichung:
dpsik = -1i / Bloch.hbar * ( Bloch.Eks(:,1) + Bloch.Eks(:,2) ) .* psi_k ...
    + 1i * hbarOmega / Bloch.hbar - Bloch.gamma / Bloch.hbar * psi_k;

% Berechnung der gesamten Polarisation:
P_t = 1 / (2 * pi)^2 * sum( conj(Bloch.dipol) .* psi_k .* Bloch.wk.' );

% Berechnung der Fourier Transformation:
fw = exp(1i * Bloch.w * t);

d_P_w = P_t * fw;
d_E_w = E_t * fw;

% Zeitentwicklung:
dpsik_E = [dpsik; d_P_w; d_E_w];
