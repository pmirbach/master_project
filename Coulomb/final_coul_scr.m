function [V , varargout] = final_coul_scr(q,para,coul_pol)

U = 1 ./ ( ( q )  .* ( 1 + para(3) .* q ) );
eps = para(2) .* ( para(2) + 1 - ( para(2) - 1 ) * exp( - para(1) .* q ) ) ./ ( para(2) + 1 + ( para(2) - 1 ) .* exp( - para(1) .* q ) );

V = U ./ eps;

V(V==inf) = coul_pol;

if nargout > 1
    varargout{1} = U;
    varargout{2} = eps;
end

return