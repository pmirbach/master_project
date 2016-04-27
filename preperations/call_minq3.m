function minq = call_minq3( Parameter, k )

b1 = Parameter.rezGV(:,1);                              % Reziproke Gittervektoren
b2 = Parameter.rezGV(:,2);

b_0 = [0; 0];                                           % Same red BZ
b_nn = [1 0 1 -1 0 -1 ; 0 1 1 0 -1 -1];                 % Next neighbours

c = [b_0, b_nn];

% All displacement vectors
b = zeros(2, size(c,2));
for ii = 1:size(c,2)
    b(1:2,ii) = c(1,ii) * b1 + c(2,ii) * b2;
end

numk = size(k,2);
numb = size(b,2);


%% Vektorisiert


allq = zeros(numk,numk,6,numb);

for nb = 1:numb
    
    [~,Kx] = meshgrid( k(1,:,1) );
    [~,Ky] = meshgrid( k(2,:,1) );
    
    Gx = repmat(b(1,nb),[numk,numk]);
    Gy = repmat(b(2,nb),[numk,numk]);
    
    for tri = 1:6
                                 
        [Kxs,~] = meshgrid( k(1,:,tri) );
        [Kys,~] = meshgrid( k(2,:,tri) );
        
        allq(:,:,tri,nb) = sqrt( ( Kx - Kxs + Gx ).^2 + ( Ky - Kys + Gy ).^2 );       
        
    end
    
end

minq = min(allq,[],4);

minq( minq < Parameter.qmin / 2 ) = 0;







