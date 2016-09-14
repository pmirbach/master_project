function lambda = get_SOC_lambda( lambda_0 , k , k_symm )

nr_k = size( k , 2 );

k_all = repmat(k,[ 1, 1, 2 ]);
k_all(:,:,1) = k_all(:,:,1) - repmat( k_symm(:,1) , [ 1, nr_k ]);
k_all(:,:,2) = k_all(:,:,2) - repmat( k_symm(:,2) , [ 1, nr_k ]);

q_all = sqrt( sum( k_all.^2 , 1 ) );
q_nearest = min( q_all , [] , 3 );

k_symm_norm = norm( k_symm(:,1) );

dummy = ( 1 - q_nearest / k_symm_norm ).^2;

lambda = lambda_0 * exp(1) * dummy .* exp( - dummy );