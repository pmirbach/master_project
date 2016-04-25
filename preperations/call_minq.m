function q_p = call_minq( Parameter, k )
    
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

% k(3,:,:) = [];

q = zeros(1,size(b,2));
q_p = zeros(size(k,2),size(k,2),size(k,3));

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
            
            q_p(nk,nks,tri) = min( q );
            
        end
        
    end
    
end


q_p( q_p < Parameter.qmin / 2 ) = 0;
