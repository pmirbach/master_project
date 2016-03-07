function [U_q] = fun_Coulomb(q,gamma,kappa)

U_q = 1 ./ ( ( q + kappa)  .* ( 1 + gamma * q ) );