function [diad_hartree, diad_fock] = fun_coul_diad(Ev_k, Ev_ks)

diad_hartree = zeros(6,6,6,6);
diad_fock = zeros(6,6,6,6);

for l1 = 1:6
    
    for l2 = 1:6
        
        diad_hartree(:,:,l1,l2) = ( conj(Ev_k(l1,:)) .* Ev_k(l1,:) ) * ...
            ( conj(Ev_ks(l2,:)) .* Ev_ks(l2,:) )';
        diad_fock(:,:,l1,l2) = ( conj(Ev_k(l1,:)) .* Ev_ks(l2,:) ) * ...
            ( conj(Ev_ks(l2,:)) .* Ev_k(l1,:) )';
        
    end
    
end