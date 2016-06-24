function CV = call_CV( Ev )
%Berechne dyadisches Produkt der Eigenvektoren, was später zur Entwicklung
%der Coulomb- und Dipolmatrixelemente in die Tight-Binding-Basis benötigt
%wird.

% kpts, tri, alpha, beta, lambda1, lambda2
CV = complex( zeros( size(Ev,3), size(Ev,4), size(Ev,1), size(Ev,1), size(Ev,2), size(Ev,2) ) );

sEv1 = size(Ev,1);
sEv2 = size(Ev,2);

for a = 1:sEv1
    
    for b = 1:sEv1
        
        for ll1 = 1:sEv2
            
            for ll2 = 1:sEv2
                
                CV( :, :, a, b, ll1, ll2) = conj( Ev( a, ll1, :, : ) ) .* Ev( b, ll2, :, : );
                
            end
            
        end
        
    end
    
end

% Tests
sample = get_sample( size( CV ) , 10 );
for ii = 1:size(sample,1)
    nk = sample(ii,1);
    ni = sample(ii,2);
    na = sample(ii,3);
    nb = sample(ii,4);
    l1 = sample(ii,5);
    l2 = sample(ii,6);
    
    CV_test = conj( Ev( na, l1, nk, ni ) ) * Ev( nb, l2, nk, ni );
    
    if ( CV( nk, ni, na, nb, l1, l2 ) ~= CV_test )
        warning('CV matrix not correct calculated!')
    end
    
end