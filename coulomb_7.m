function coulomb_7(constAg,Para,Data,Prep)


for nll = 1:size(Para.coul_indices,1)
    
    l1 = Para.coul_indices(nll,1);
    l2 = Para.coul_indices(nll,2);
    
    if l1 <= 3
        d = 0;
        
        e = 0;
        f = 3;
        h = 3;
    else
        d = 3;
        
        e = 3;
        f = 0;
        h = -3;
    end
    
    fprintf('l1 = %d, l2 = %d\n',[l1,l2+h])
    
    for a = 1:3
        
        for b = 1:3
            
                
%                 fprintf('l1 = %d, l2 = %d, a = %d, b = %d\n',[l1,l2+h,a+e,b+f])
                
%                 V_hartree_off(:,:,nll) = V_hartree_off(:,:,nll) + ...
%                     ( Prep.CV(:,1,a+e,a+e,l1,l1) * Prep.CV(:,tri,b+f,b+f,l2+h,l2+h).' ) * V_orbital_h(a,b);
                               
            
        end
        
    end
    
end
