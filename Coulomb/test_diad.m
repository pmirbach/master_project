function test_diad(Ev_k, Ev_ks,combi)

% Fock artig
% combi = [1 1 1 1];

diad_einzeln = zeros(6);

for a = 1:6
    
    for b = 1:6
        
        diad_einzeln(a,b) = conj(Ev_k(a,combi(1))) * conj(Ev_ks(b,combi(2))) *...
            Ev_k(a,combi(3)) * Ev_ks(b,combi(4));
        
    end
    
end

disp(diad_einzeln)

% diad_cool = (conj(Ev_k(:,1)) * (conj(Ev_ks(:,1)))') .* (Ev_k(:,1) * (Ev_ks(:,1))');