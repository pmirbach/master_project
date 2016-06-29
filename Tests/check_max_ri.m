function check_max_ri( A )

Areal = abs( real( A ) );
Aimag = abs( imag( A ) );

fprintf( 'Maximaler Realteil    : %e \n', max( Areal(:) ) )
fprintf( 'Maximaler Imaginärteil: %e \n', max( Aimag(:) ) )