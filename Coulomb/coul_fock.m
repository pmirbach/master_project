function [coul_diad_f] = coul_fock(Ev_k, Ev_ks, l1, l2)

coul_diad_f = zeros(6,6,6);

for a = 1:6
    
    for b = 1:6
        
        for ni = 1:6
            
            coul_diad_f(a,b,ni) = conj(Ev_k(a,l1)) *...
                conj(Ev_ks(b,l2,ni)) * ...
                Ev_k(b,l1) * Ev_ks(a,l2,ni) ;
        end
        
    end
    
end

end