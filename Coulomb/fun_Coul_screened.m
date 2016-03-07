function [V] = fun_Coul_screened(q,para,kappa)

d = para(1);
eps_inf = para(2);
gamma = para(3);

diel = eps_inf * ( eps_inf + 1 - (eps_inf - 1) * exp(-d * q) ) ...
    ./ ( eps_inf + 1 + (eps_inf - 1) * exp(-d * q) );

U = 1 ./ ( ( q + kappa)  .* ( 1 + gamma * q ) );

V = U ./ diel;