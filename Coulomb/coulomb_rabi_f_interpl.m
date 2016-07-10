function [V_rabi_fock] = coulomb_rabi_f_interpl(Ctrl, Para, Prep )

V_rabi_fock = zeros( Para.nr.k, Para.nr.k, size(Para.dipol_trans,1) );

para_map = [1 2 3 ; 2 4 5 ; 3 5 6];

for nll = 1:Para.nr.dipol
    
    hh = Para.dipol_trans(nll,1);
    ee = Para.dipol_trans(nll,2);
    
    if hh <= 3
        d = 0;
    else
        d = 3;
    end
    
    for tri = 1:6
        
        Coul = 0;
        
        for a = 1:3
            
            for b = 1:3
                
%                 reshape( Prep.V_ab_interpl{ para_map(a,b) } ( reshape( Prep.minq(:,:,tri) , 1 , [] ) ), size( Prep.minq(:,:,tri) ) )
                             
                Coul = Coul + ( Prep.CV_noSOC(:,1,a,b,ee-d,hh-d) * Prep.CV_noSOC(:,tri,b,a,hh-d,ee-d).' ) .* ...
                    reshape( Prep.V_ab_interpl{ para_map(a,b) } ( reshape( Prep.minq(:,:,tri) , 1 , [] ) ), size( Prep.minq(:,:,tri) ) );

            end
            
        end
        
        V_rabi_fock(:,:,nll) = V_rabi_fock(:,:,nll) + abs( Coul );
        
    end
    
end


V_rabi_fock = Para.vorf.coul * abs( V_rabi_fock );
