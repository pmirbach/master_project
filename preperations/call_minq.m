function minq = call_minq( Para, k )

allq = zeros( Para.nr.k , Para.nr.k , 6 , Para.nr.b );

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

% Tests
sample = get_sample( size( minq ) , 10 );
for ii = 1:size(sample,1)
    nk = sample(ii,1);
    nks = sample(ii,2);
    tri = sample(ii,3);
    
    minq_test = zeros(1,Para.nr.b);
    for nb = 1:Para.nr.b
        minq_test(nb) = norm( k(:,nk,1) - k(:,nks,tri) + Para.k.b(:,nb) );
    end
    
    minq_test( minq_test < 100 * eps ) = 0;
    
    if ( min( minq_test ) - minq(nk,nks,tri) > 100 * eps )
        warning('Minimal q matrix not correct calculated!')
    end
    
end
