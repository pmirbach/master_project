function b = call_bnn( GV )

b1 = GV(:,1);                              % Reziproke Gittervektoren
b2 = GV(:,2);

b_0 = [0; 0];                                           % Same red BZ
b_nn = [1 0 1 -1 0 -1 ; 0 1 1 0 -1 -1];                 % Next neighbours

c = [b_0, b_nn];

% All displacement vectors
b = zeros(2, size(c,2));
for ii = 1:size(c,2)
    b(1:2,ii) = c(1,ii) * b1 + c(2,ii) * b2;
end