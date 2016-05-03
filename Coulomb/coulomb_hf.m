function [V_fock, V_hartree, V_hartree_off] = coulomb_hf(Ctrl, Para, Prep)

V_fock = zeros( Para.nr.k, Para.nr.k, size(Para.coul_indices,1) );
V_hartree = V_fock;
V_hartree_off = V_fock;


if Ctrl.carrier_density > 0
    
    V_orbital_h = final_coul_long(Para.coul.screened);
    
    para_map = [1 2 3 ; 2 4 5 ; 3 5 6];
    para_map = repmat(para_map,[2,2]);
    
    for nll = 1:size(Para.coul_indices,1)
        
        l1 = Para.coul_indices(nll,1);
        l2 = Para.coul_indices(nll,2);
        l3 = Para.coul_indices(nll,3);
        
        if l1 <= 3
            d = 0;
            
            e = 0;
            f = 3;
        else
            d = 3;
            
            e = 3;
            f = 0;
        end
        
        for a = 1:3
            
            for b = 1:3
                
                for tri = 1:6
                    
                    V_fock(:,:,nll) = V_fock(:,:,nll) + ...
                        ( Prep.CV(:,1,a+d,b+d,l1,l1) * Prep.CV(:,tri,b+d,a+d,l2,l2).' ) .* ...
                        final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                    
                    V_hartree(:,:,nll) = V_hartree(:,:,nll) + ...
                        ( Prep.CV(:,1,a+d,a+d,l1,l1) * Prep.CV(:,tri,b+d,b+d,l2,l2).' ) * V_orbital_h(a,b);
                    
                    V_hartree_off(:,:,nll) = V_hartree_off(:,:,nll) + ...
                        ( Prep.CV(:,1,a+e,a+e,l1,l1) * Prep.CV(:,tri,b+f,b+f,l3,l3).' ) * V_orbital_h(a,b);
                    
                end
                
            end
            
        end
        
    end
    V_fock = Para.vorf.coul * abs( V_fock );
    V_hartree = Para.vorf.coul * ( V_hartree );
    V_hartree_off = Para.vorf.coul * ( V_hartree_off );
    
end
