function [V_rabi_fock] = coulomb_rabi_f_2(Ctrl, Para, Prep)


V_rabi_fock = zeros( Para.nr.k, Para.nr.k, Para.nr.tri, size(Para.dipol_trans,1) );

para_map = [1 2 3 ; 2 4 5 ; 3 5 6];

for nll = 1:Para.nr.dipol
    
    hh = Para.dipol_trans(nll,1);
    ee = Para.dipol_trans(nll,2);
    
    if hh <= 3
        d = 0;
    else
        d = 3;
    end
    
    for a = 1:3
        
        for b = 1:3
            
            for tri = 1:6
                      
                V_rabi_fock(:,:,tri,nll) = V_rabi_fock(:,:,tri,nll) + ...
                    ( Prep.CV_noSOC(:,1,a,b,ee-d,hh-d) * Prep.CV_noSOC(:,tri,b,a,hh-d,ee-d).' ) .* ...
                    final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                                   
            end
            
        end
        
    end
    
end

V_rabi_fock = abs( V_rabi_fock );
