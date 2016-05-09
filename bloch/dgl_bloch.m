function dpsik_E = dgl_bloch(t,psik_E,Bloch)

% Nur Übergang 1 -> 2

psi_k = psik_E(1:Bloch.nrk);


E_t = Bloch.E0 * exp(-1/2 * ( ( t - Bloch.t_peak ) / Bloch.sigma )^2 * 4 * log(2) );


dpsik = -1i / Bloch.hbar * ( Bloch.Eks(1,:) + Bloch.Eks(2,:) )' .* psi_k ...
    + 1i * E_t * conj(Bloch.dipol) / Bloch.hbar - Bloch.gamma / Bloch.hbar * psi_k;

% dpsik = - gamma / constAg.hbar * psi_k;

% P_t = 1 / (2 * pi)^2 * sum( conj(Bloch.dipol) .* psi_k );
P_t = 1 / (2 * pi)^2 * sum( conj(Bloch.dipol) .* psi_k .* Bloch.wk.' );



fw = exp(1i * Bloch.w * t);

d_P_w = P_t * fw;
d_E_w = E_t * fw;


dpsik_E = [dpsik; d_P_w; d_E_w];
