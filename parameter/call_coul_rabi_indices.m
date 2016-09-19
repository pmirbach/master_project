function [ coul_rabi_unique , coul_dip_mapping ] = call_coul_rabi_indices( dipol_trans )

unique_coul = dipol_trans;

unique_coul( all( unique_coul > 3, 2 ) , : ) = unique_coul( all( unique_coul > 3, 2 ) , : ) - 3;


[coul_rabi_unique, ~, coul_dip_mapping] = unique(unique_coul,'rows');
