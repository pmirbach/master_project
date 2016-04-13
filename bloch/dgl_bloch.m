function dpsik_E = dgl_bloch(t,psik_E,Data,constAg)

% Nur Übergang 1 -> 2

psi_k = psik_E(1:631);

E0 = 1e-6;
t_peak = 0.003;
sigma = 0.001;


E_t = E0 * exp(-1/2 * ( ( t - t_peak ) / sigma )^2 * 4 * log(2) );

d = 1 / sqrt(2) * transpose(Data.dipol{2,1}(1,:) - 1i * Data.dipol{2,1}(2,:));

gamma = 10;

dpsik = -1i / constAg.hbar * ( Data.Ek_s(1,:) + Data.Ek_s(2,:) )' .* psi_k ...
    + 1i * E_t * conj(d) / constAg.hbar - gamma / constAg.hbar * psi_k;

% dpsik = - gamma / constAg.hbar * psi_k;

P_t = 1 / (2 * pi)^2 * sum( conj(d) .* psi_k );

w = (-1000:1000)';
fw = exp(1i * w * t);

d_P_w = P_t * fw;
d_E_w = E_t * fw;


dpsik_E = [dpsik; d_P_w; d_E_w];
