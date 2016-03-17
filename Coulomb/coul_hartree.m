function [coul_diad_h] = coul_hartree(Ev_k, Ev_ks, l1, l2)

coul_diad_h = zeros(6,6,6);

for a = 1:6
    
    for b = 1:6
        
        for ni = 1:6
            
            coul_diad_h(a,b,ni) = conj(Ev_k(a,l1)) * ...
                conj(Ev_ks(b,l2,ni)) * ...
                Ev_ks(b,l2,ni) * Ev_k(a,l1) ;
           
        end
        
    end
    
end

end