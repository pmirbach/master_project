function [ V_rabi_fock ] = coulomb_rabi_f_interpl(Ctrl, Para, Prep )

if Ctrl.Coul.active
    
    V_rabi_fock = zeros( Para.nr.k, Para.nr.k, Para.nr.dipol_coul );
    
    para_map = [1 2 3 ; 2 4 5 ; 3 5 6];
    
    for nll = 1:Para.nr.dipol_coul
        
        hh = Para.coul_rabi_unique(nll,1);
        ee = Para.coul_rabi_unique(nll,2);
                
        for tri = 1:6
            
            Coul = 0;
            
            for a = 1:3
                
                for b = 1:3
                    
                    Coul = Coul + ( Prep.CV_noSOC(:,1,a,b,ee,hh) * Prep.CV_noSOC(:,tri,b,a,hh,ee).' ) .* ...
                        reshape( Prep.V_ab_interpl{ para_map(a,b) } ( reshape( Prep.minq(:,:,tri) , 1 , [] ) ), size( Prep.minq(:,:,tri) ) );
                    
                end
                
            end
            
            V_rabi_fock(:,:,nll) = V_rabi_fock(:,:,nll) + abs( Coul );
            
        end
        
    end
    
    V_rabi_fock = Para.vorf.coul * abs( V_rabi_fock );
    
else
    V_rabi_fock = [];
end
