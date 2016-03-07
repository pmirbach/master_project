function [diel] = fun_diel2(q, d, eps )

diel = eps * ( eps + 1 - (eps - 1) * exp(-d * q) ) ...
    ./ ( eps + 1 + (eps - 1) * exp(-d * q) );