function [eps_q] = fun_diel(q, eps, d)


eps_q = eps * ( eps + 1 - (eps - 1) .* exp(-d * q) ) ...
    ./ ( eps + 1 + (eps - 1) .* exp(-d * q) );
