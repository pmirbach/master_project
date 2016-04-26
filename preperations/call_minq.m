function minq = call_minq( Parameter, k )

b1 = Parameter.rezGV(:,1);                              % Reziproke Gittervektoren
b2 = Parameter.rezGV(:,2);

b_0 = [0; 0];                                           % Same red BZ
b_nn = [1 0 1 -1 0 -1 ; 0 1 1 0 -1 -1];                 % Next neighbours

c = [b_0, b_nn];

% All displacement vectors
b = zeros(2, size(c,2));
for nb = 1:size(c,2)
    b(:,nb) = c(1,nb) * b1 + c(2,nb) * b2;
end

numk = size(k,2);
numb = size(b,2);

% k(3,:,:) = [];

q = zeros(1,size(b,2));
minq = zeros(size(k,2),size(k,2),size(k,3));

for nk = 1:numk
    
    %     disp(nk)
    
    kx = k(1,nk,1);
    ky = k(2,nk,1);
    
    for nks = 1:numk
        
        for tri = 1:6
            
            kxs = k(1,nks,tri);
            kys = k(2,nks,tri);
            
            for nb = 1:numb
                
                q(nb) = sqrt( (kx - kxs + b(1, nb) )^2 + (ky - kys + b(2, nb))^2 );
                
            end
            
            minq(nk,nks,tri) = min( q );
            
        end
        
    end
    
end


minq( minq < Parameter.qmin / 2 ) = 0;





% minq2 = 1e5 * ones(numk,numk,6);

allq = zeros(numk,numk,6,numb);

[Kx,~] = meshgrid(k(1,:,1));
[Ky,~] = meshgrid(k(2,:,1));

for tri = 3
    
    [~,Kxs] = meshgrid(k(1,:,tri));
    [~,Kys] = meshgrid(k(2,:,tri));
    
    
    for nb = 1:numb
        
        Gx = repmat(b(1,nb),[numk,numk]);
        Gy = repmat(b(2,nb),[numk,numk]);
        
        allq(:,:,tri,nb) = sqrt( ( Kx - Kxs + Gx ).^2 + ( Ky - Kys + Gy ).^2 );
        
%                 1
    end
    
end

minq2 = min(allq,[],4);
minq3 = permute(minq2,[2,1,3]);

max(max(max(abs(minq3-minq))))

1








