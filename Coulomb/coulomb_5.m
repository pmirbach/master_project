function [V_fock, V_hartree] = coulomb_5(constAg,Para,Data,Prep)



% % % % % % % % %  V orbital h hier berechnen
% fun_coul_orbital_hartree <----

V_orbital_h = fun_coul_orbital_hartree(Para.coul.screened);

para_map = [1 2 3 ; 2 4 5 ; 3 5 6];
para_map = repmat(para_map,[2,2]);


V_fock = zeros( Para.nr.k, Para.nr.k, size(Para.coul_indices,1) );
V_hartree = V_fock;

% tic

for nll = 1:size(Para.coul_indices,1)
    
    l1 = Para.coul_indices(nll,1);
    l2 = Para.coul_indices(nll,2);
    
    if l1 <= 3
        d = 0;
    else
        d = 3;
    end
    
    for a = 1:3
        
        for b = 1:3
            
            for tri = 1:6
                                
                V_fock(:,:,nll) = V_fock(:,:,nll) + ...
                    abs( Prep.CV(:,1,a+d,b+d,l1,l1) * Prep.CV(:,tri,b+d,a+d,l2,l2).' ) .* ...
                    final_coul_scr(Prep.minq(:,:,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                
                V_hartree(:,:,nll) = V_hartree(:,:,nll) + ...
                    abs( Prep.CV(:,1,a+d,a+d,l1,l1) * Prep.CV(:,tri,b+d,b+d,l2,l2).' ) * V_orbital_h(a,b);
                
                
            end
            
        end
        
    end
    
end
V_fock = Para.vorf.coul * V_fock;
V_hartree = Para.vorf.coul * V_hartree;

% toc

% 1
% 
% figure
% hold on
% % for ii = 1:6
% scatter3(Data.k(1,:,1),Data.k(2,:,1),V_fock(1,:,1)')
% % end
% 
% 1
