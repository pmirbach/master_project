function [kneu] = umklapp2(Parameter,k,k0)

b1 = Parameter.rezGV(:,1);          % Reziproke Gittervektoren
b2 = Parameter.rezGV(:,2);

b_0 = [0; 0];                       % Same red BZ
b_nn = [1 0 1 -1 0 -1 ; ...         % Next neighbours
    0 1 1 0 -1 -1];

c = [b_0, b_nn];

kneu = zeros(size(k));

% All displacement vectors
b = zeros(2, size(c,2));            
for ii = 1:size(c,2)
    b(1:2,ii) = c(1,ii) * b1 + c(2,ii) * b2;
end

for nk = 1:size(k,2)
    
    for ntri = 1:size(k,3)
        
        [M,I] = min(sqrt(sum( (repmat( k(:,nk,ntri) - k0 ,1,size(b,2)) + b).^2,1)));
        kneu(:,nk,ntri) = k(:,nk,ntri) + b(:,I);
        
    end
    
end











