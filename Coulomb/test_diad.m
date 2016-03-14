function [diad_h, diad_f] = test_diad(Ev_k, Ev_ks,combi)

% Fock artig
% combi = [1 1 1 1];

diad_h = zeros(6);
diad_f = zeros(6);

for a = 1:6
    
    for b = 1:6
        
        diad_h(a,b) = conj(Ev_k(a,combi(1))) * conj(Ev_ks(b,combi(2))) *...
            Ev_ks(b,combi(3)) * Ev_k(a,combi(4));
        
        diad_f(a,b) = conj(Ev_k(a,combi(1))) * conj(Ev_ks(b,combi(2))) *...
            Ev_k(b,combi(3)) * Ev_ks(a,combi(4));
        
    end
    
end

% disp(diad_einzeln)

% diad_cool = (conj(Ev_k(:,1)) * (conj(Ev_ks(:,1)))') .* (Ev_k(:,1) * (Ev_ks(:,1))');