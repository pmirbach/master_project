function [V_orbital_f] = fun_coul_orbital_fock2(q,para,kappa)

V_temp = zeros(6);
V_orbital_f = zeros(3,3,6);

for jj = 1:size(para,1)
   
%     d = para(jj,1);
%     eps_inf = para(jj,2);
%     gamma = para(jj,3);
    
    diel = para(jj,2) * ( para(jj,2) + 1 - (para(jj,2) - 1) * exp(-para(jj,1) * q) ) ...
        ./ ( para(jj,2) + 1 + (para(jj,2) - 1) * exp(-para(jj,1) * q) );
    
    V_temp(jj,:) = 1 ./ ( ( q + kappa)  .* ( 1 + para(jj,3) * q ) ) ./ diel;
    
end

V_orbital_f(1,1,:) = V_temp(1,:);
V_orbital_f(1,2,:) = V_temp(2,:);
V_orbital_f(2,1,:) = V_temp(2,:);
V_orbital_f(1,3,:) = V_temp(3,:);
V_orbital_f(3,1,:) = V_temp(3,:);
V_orbital_f(2,2,:) = V_temp(4,:);
V_orbital_f(2,3,:) = V_temp(5,:);
V_orbital_f(3,2,:) = V_temp(5,:);
V_orbital_f(3,3,:) = V_temp(6,:);



V_orbital_f = repmat(V_orbital_f,2,2);