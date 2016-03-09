function [V_hartree, V_fock] = fun_coul_orbital_m(q,para,kappa)


V_hartree_s = zeros(1,6);
V_fock_s = zeros(1,6);

for jj = 1:size(para,1)
    
    V_fock_s(jj) = sum( fun_Coul_screened(q,para(jj,:),kappa) );
    
    V_hartree_s(jj) = 6 * fun_Coul_screened_long(para(jj,:));
    
end

V_hartree_spin = [V_hartree_s(1), V_hartree_s(2), V_hartree_s(3); 
    V_hartree_s(2), V_hartree_s(4), V_hartree_s(5); 
    V_hartree_s(3), V_hartree_s(5), V_hartree_s(6)];

V_fock_spin = [V_fock_s(1), V_fock_s(2), V_fock_s(3); 
    V_fock_s(2), V_fock_s(4), V_fock_s(5); 
    V_fock_s(3), V_fock_s(5), V_fock_s(6)];

V_hartree = blkdiag(V_hartree_spin,V_hartree_spin);
V_fock = blkdiag(V_fock_spin, V_fock_spin);