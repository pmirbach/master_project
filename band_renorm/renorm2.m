function [Ek_hf, Ek_h, Ek_f] = renorm2(Para, Ek, V_f, V_h, fk, wk, ll)


Ek_h = Ek(:,:,1);
Ek_f = Ek(:,:,1);
Ek_hf = Ek(:,:,1);

gew = wk * Para.BZsmall.area / 6;




coul_map = reshape(1:9,[3,3])';
coul_map = blkdiag(coul_map,coul_map+9);




vorz = [-1 1 1; 1 -1 -1; 1 -1 -1];
vorz = blkdiag(vorz,vorz);

for nll = 1:size(ll,1)
    
    l1 = ll(nll,1);
    l2 = ll(nll,2);
   
    Ek_f(l1,:) = Ek_f(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( fk(l2,:) .* gew ) * (-1) * V_f(:,:,coul_map(l1,l2)).';
%     Ek_h(l1,:) = Ek_h(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( fk(l2,:) .* wk / 6 ) * V_h(:,:,l2)';
%     Ek_hf(l1,:) = Ek_hf(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( fk(l2,:) .* wk / 6 ) * ( V_h(:,:,l2) - V_f(:,:,l2))';
    
end

