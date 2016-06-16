function [nEv] = EV_ortho( Para, k , Ek , Ev )

nEv = 1;

kn1 = repmat( k , 3 , 1 );
kn = reshape( kn1 , 2, [] );

Ekn = reshape( Ek, 1, [] );
Ekv = reshape( Ev, 3, [] );


[B,I] = sort(Ekn);


for ii = 1:size(B,2) - 1
    diff = B(ii+1) - B(ii);
    if diff < 1e-16
        
        ind1 = I(ii);
        ind2 = I(ii+1);
        
%         Ekv(:,ind2)' * Ekv(:,ind1)
        
        k1 = kn(:,ind1);
        k2 = kn(:,ind2);
        
        H_TB1 = TB_Liu_TNN_fun(k1, Para.TB);
        H_TB2 = TB_Liu_TNN_fun(k2, Para.TB);
        
        
        H_TB1*1e3 * Ekv(:,ind1) - Ekn(ind1) * Ekv(:,ind1)
        H_TB2*1e3 * Ekv(:,ind2) - Ekn(ind2) * Ekv(:,ind2)
        
        
        
        OrEv = orth( [ Ekv(:,ind1) , Ekv(:,ind2) ] );
        
       
        v2 = Ekv(:,ind2) - Ekv(:,ind1)' * Ekv(:,ind2) / ( Ekv(:,ind1)' * Ekv(:,ind1) ) * Ekv(:,ind1);
%         OrEv2 = [Ekv(:,ind1) , v2];
        
        H_TB1*1e3 * OrEv(:,1) - Ekn(ind1) * OrEv(:,1)
        H_TB2*1e3 * OrEv(:,2) - Ekn(ind2) * OrEv(:,2)   
        
        H_TB2*1e3 * v2 - Ekn(ind2) * v2  
        
        1
        
    end
end


1

% Ekn2 = reshape( Ekn , size(Ek) );




