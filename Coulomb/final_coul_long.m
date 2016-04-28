function v_long_matrix = final_coul_long(para)

nr_orbitals = size(para,1);

v_long_matrix = zeros(3);
v_temp = zeros(1,6);

for jj = 1:nr_orbitals

    v_temp(jj) = -( para(jj,3) + para(jj,4) );
    
end


v_long_matrix(1,1) = v_temp(1);
v_long_matrix(1,2) = v_temp(2);
v_long_matrix(2,1) = v_temp(2);
v_long_matrix(1,3) = v_temp(3);
v_long_matrix(3,1) = v_temp(3);
v_long_matrix(2,2) = v_temp(4);
v_long_matrix(2,3) = v_temp(5);
v_long_matrix(3,2) = v_temp(5);
v_long_matrix(3,3) = v_temp(6);
