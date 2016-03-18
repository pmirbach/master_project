function [V_orbital_f] = fun_coul_orbital_fock(q,para,kappa)

ind = {[1,1],[1,2],[1,3],[2,2],[2,3],[3,3]};

V_orbital_f = zeros(3,3,6);

for jj = 1:size(para,1)
    
    V_orbital_f(ind{jj}(1),ind{jj}(2),:) = ...
        fun_Coul_screened(q,para(jj,:),kappa);
    
    V_orbital_f(ind{jj}(2),ind{jj}(1),:) = ...
        fun_Coul_screened(q,para(jj,:),kappa);
    
end

V_orbital_f = repmat(V_orbital_f,2,2);