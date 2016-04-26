function [V_fock, V_hartree] = coulomb_5(constAg,Parameter,Data,Prep)

numk = size(Data.k,2);

vorf = constAg.ec^2 / ( 2 * constAg.eps_0);


[A,B] = meshgrid(1:3,1:3);
c=cat(2,A,B);
ll=reshape(c,[],2);
ll = [ll; ll+3];

para_map = [1 2 3 ; 2 4 5 ; 3 5 6];


V_fock = zeros( numk, numk, size(ll,1) );
V_hartree = V_fock;

d = 0;
% tic

for nll = 1:size(ll,1)
    
    l1 = ll(nll,1);
    l2 = ll(nll,2);
    
%     if l1 == 4
%         d = 3;
%     end
    
    for a = 1:3
        
        for b = 1:3
            
            for tri = 1:6
                                
                V_fock(:,:,nll) = V_fock(:,:,nll) + vorf * ...
                    abs( Prep.CV2(:,1,a+d,b+d,l1,l1) * Prep.CV2(:,tri,b+d,a+d,l2,l2).' ) .* ...
                    final_coul_scr(Prep.minq(:,:,tri),Parameter.coul_screened(para_map(a,b),:),Parameter.coul_pol);
                
                V_hartree(:,:,nll) = V_hartree(:,:,nll) + vorf * ...
                    abs( Prep.CV2(:,1,a+d,a+d,l1,l1) * Prep.CV2(:,tri,b+d,b+d,l2,l2).' ) * Prep.V_orbital_h(a,b);
                
            end
            
        end
        
    end
    
    
end

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
