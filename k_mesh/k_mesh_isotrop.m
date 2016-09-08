function [ k, w_k , n_k , index_iso ] = k_mesh_isotrop( q_max, N_q, N_phi )


q = linspace(0,q_max,N_q+1);
q(1) = [];

phi = 0:2*pi/(N_phi):2*pi-0.0001;

index_iso = 2:N_phi:(N_q*N_phi+1);
index_iso = [index_iso; index_iso + N_phi - 1];

k = zeros(2 , N_phi * N_q + 1);

for ii = 1:N_q
    k(1,index_iso(1,ii):index_iso(2,ii)) = cos(phi) * q(ii);
    k(2,index_iso(1,ii):index_iso(2,ii)) = sin(phi) * q(ii);
end

n_k = size(k,2);

w_k = 1;