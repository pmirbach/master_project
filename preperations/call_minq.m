function minq = call_minq( Para, k )

allq = zeros(Para.nr.k,Para.nr.k,6,Para.nr.b);

for nb = 1:Para.nr.b
    
    [~,Kx] = meshgrid( k(1,:,1) );
    [~,Ky] = meshgrid( k(2,:,1) );
    
    Gx = repmat(Para.k.b(1,nb),[Para.nr.k,Para.nr.k]);
    Gy = repmat(Para.k.b(2,nb),[Para.nr.k,Para.nr.k]);
    
    for tri = 1:6
                                 
        [Kxs,~] = meshgrid( k(1,:,tri) );
        [Kys,~] = meshgrid( k(2,:,tri) );
        
        allq(:,:,tri,nb) = sqrt( ( Kx - Kxs + Gx ).^2 + ( Ky - Kys + Gy ).^2 );       
        
    end
    
end

minq = min(allq,[],4);

minq( minq < Para.k.qmin / 2 ) = 0;

