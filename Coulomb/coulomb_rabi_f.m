function [V_rabi_fock] = coulomb_rabi_f(Ctrl, Para, Prep)

V_rabi_fock = zeros( Para.nr.k, Para.nr.k, size(Para.dipol_trans,1) );

para_map = [1 2 3 ; 2 4 5 ; 3 5 6];

for nll = 1:size(Para.dipol_trans,1)
    
    l1 = Para.dipol_trans(nll,1);
    l2 = Para.dipol_trans(nll,2);
    
    if l1 <= 3
        d = 0;
    else
        d = 3;
    end
    
    for a = 1:3
        
        for b = 1:3
            
            for tri = 1:6
                
                V_rabi_fock(:,:,nll) = V_rabi_fock(:,:,nll) + ...
                    ( Prep.CV(:,1,a+d,b+d,l2,l1) * Prep.CV(:,tri,b+d,a+d,l1,l2).' ) .* ...
                    final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                
            end
            
        end
        
    end
    
end

V_rabi_fock = Para.vorf.coul * abs( V_rabi_fock );