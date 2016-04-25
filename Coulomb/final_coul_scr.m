function V = final_coul_scr(q,para,coul_pol)

V = 1 ./ ( ( q )  .* ( 1 + para(3) .* q ) ) ./ ...
    para(2) .* ( para(2) + 1 - (para(2) - 1) * exp(-para(1) .* q) ) ./ ( para(2) + 1 + (para(2) - 1) .* exp(-para(1) .* q) );

V(V==inf) = coul_pol;

return