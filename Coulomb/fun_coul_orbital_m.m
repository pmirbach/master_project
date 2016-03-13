function [V_orbital_h, V_orbital_f] = fun_coul_orbital_m(q,para,kappa)

ind_l = [1 2 3 5 6 9];

V_orbital_h = zeros(3);
V_orbital_f = zeros(3);

for jj = 1:size(para,1)
    V_orbital_h(ind_l(jj)) = 6 * fun_Coul_screened_long(para(jj,:));
    V_orbital_f(ind_l(jj)) = sum( fun_Coul_screened(q,para(jj,:),kappa) );
end

V_orbital_h = V_orbital_h + V_orbital_h' - diag(diag(V_orbital_h));
V_orbital_f = V_orbital_f + V_orbital_f' - diag(diag(V_orbital_f));

V_orbital_h = blkdiag(V_orbital_h, V_orbital_h);
V_orbital_f = blkdiag(V_orbital_f, V_orbital_f);