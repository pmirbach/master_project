function [V_orbital_h] = fun_coul_orbital_hartree(para)

ind = {1,[2,4],[3,7],5,[6,8],9};

V_orbital_h = zeros(3);

for jj = 1:size(para,1)
    
    V_orbital_h(ind{jj}) = fun_Coul_screened_long(para(jj,:));
    
end

V_orbital_h = repmat(V_orbital_h,2,2);