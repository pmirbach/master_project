function [coul_intrp] = call_coul_intrp(constAg, Parameter)

q = 0.0:0.0001:14;

vorf = constAg.ec^2 / ( 2 * constAg.eps_0);

coul_intrp_tmp = cell(1,6);

for jj = 1:6
    
    para = Parameter.coul_screened(jj,:);

    [V,U] = fun_Coul_screened(q, para, Parameter.coul_kappa);
    
    V = vorf * V;
        
    coul_intrp_tmp{jj} = griddedInterpolant(q,V);

end

coul_intrp = cell(3);
coul_intrp{1,1} = coul_intrp_tmp{1};
coul_intrp{1,2} = coul_intrp_tmp{2};
coul_intrp{2,1} = coul_intrp_tmp{2};
coul_intrp{1,3} = coul_intrp_tmp{3};
coul_intrp{3,1} = coul_intrp_tmp{3};
coul_intrp{2,2} = coul_intrp_tmp{4};
coul_intrp{2,3} = coul_intrp_tmp{5};
coul_intrp{3,2} = coul_intrp_tmp{5};
coul_intrp{3,3} = coul_intrp_tmp{6};

