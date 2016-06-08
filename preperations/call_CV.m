function CV = call_CV( Ev )
%Berechne dyadisches Produkt der Eigenvektoren, was später zur Entwicklung
%der Coulomb- und Dipolmatrixelemente in die Tight-Binding-Basis benötigt
%wird.

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

% c Daniel

end