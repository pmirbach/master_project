function [V_fock, V_hartree] = coulomb_noSOC(Para,Prep)

V_long = final_coul_long(Para.coul.screened);

para_map = [1 2 3 ; 2 4 5 ; 3 5 6];

V_fock = zeros( Para.nr.k, Para.nr.k, size(Para.coul_indices,1) / 2 );
V_hartree = V_fock;


for nll = 1:size(Para.coul_indices,1) / 2
    
    l1 = Para.coul_indices(nll,1);
    l2 = Para.coul_indices(nll,2);
        
    for a = 1:3
        
        for b = 1:3
            
            for tri = 1:6
                                
                V_fock(:,:,nll) = V_fock(:,:,nll) + ...
                    ( Prep.CV_noSOC(:,1,a,b,l1,l1) * Prep.CV_noSOC(:,tri,b,a,l2,l2).' ) .* ...
                    final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                
                V_hartree(:,:,nll) = V_hartree(:,:,nll) + ...
                    ( Prep.CV_noSOC(:,1,a,a,l1,l1) * Prep.CV_noSOC(:,tri,b,b,l2,l2).' ) * V_long(a,b);
                
                
            end
            
        end
        
    end
    
end

V_fock = Para.vorf.coul * abs( V_fock );
V_hartree = Para.vorf.coul * ( V_hartree );

