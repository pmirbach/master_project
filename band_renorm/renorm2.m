function [Ek_hf, Ek_h, Ek_f, ren_h] = renorm2(Para, Ek, V_f, V_h, V_h_off, fk, wk, ll)


Ek_h = Ek(:,:,1);
Ek_f = Ek(:,:,1);
Ek_hf = Ek(:,:,1);

gew = wk * Para.BZsmall.area / 6;




coul_map = reshape(1:9,[3,3])';
coul_map = repmat(coul_map,[2,2]);
coul_map(4:6,:) = coul_map(4:6,:) + 9;
% coul_map = blkdiag(coul_map,coul_map+9);
% coul_map = blkdiag(coul_map,coul_map);

ren_h = zeros(size(Ek_h));

vorz = [-1 1 1; -1 1 1; -1 1 1];
vorz = blkdiag(vorz,vorz);

for nll = 1:size(ll,1)
    
    l1 = ll(nll,1);
    l2 = ll(nll,2);
    l3 = ll(nll,3);
    
%     ren_hartree = ( ( fk(l2,:) + fk(l3,:) ) .* gew ) * V_h(:,:,coul_map(l1,l2))';
    
    ren_hartree = ( fk(l2,:) .* gew ) * V_h(:,:,coul_map(l1,l2)).';
    
%     ren_hartree = ( fk(l2,:) .* gew ) * V_h(:,:,coul_map(l1,l2)).' + ( fk(l3,:) .* gew ) * V_h_off(:,:,coul_map(l1,l3)).' ;
%     ren_hartree = ( fk(l3,:) .* gew ) * V_h_off(:,:,coul_map(l1,l2)).' ;

    ren_h(l1,:) = ren_h(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ren_hartree;
 
    
    
    
    ren_fock = - ( fk(l2,:) .* gew ) * V_f(:,:,coul_map(l1,l2)).';
    
    Ek_f(l1,:) = Ek_f(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ren_fock;

    Ek_h(l1,:) = Ek_h(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ren_hartree;
    
    Ek_hf(l1,:) = Ek_hf(l1,:) + 1 / ( 2 * pi )^2 * vorz(l1,l2) * ( ren_hartree + ren_fock );
    
end

 