function b = call_bnn( GV )

b1 = GV(:,1);                              % Reziproke Gittervektoren
b2 = GV(:,2);

b_nn = [ [ 0 ; 0 ] , [ 0 ; 1 ] , [ 0 ; -1 ] , ...
    [ 1 ; 0 ] ,  [ 1 ; 1 ] , [ 1 ; -1 ] , ...
    [ -1 ; 0 ],  [ -1 ; 1 ] , [ -1 ; -1 ] ];       % Rec lattice vectors until next neighbours

% All displacement vectors
b_tmp = zeros( 2 , size( b_nn , 2 ) );

for ii = 1:size( b_nn , 2 )
    b_tmp(1:2,ii) = b_nn(1,ii) * b1 + b_nn(2,ii) * b2;
end

b_tmp_norm = sqrt( sum( b_tmp .^2 , 1 ) );

[ ~, b_ind ] = sort( b_tmp_norm );

b = b_tmp( : , b_ind(1:7) );