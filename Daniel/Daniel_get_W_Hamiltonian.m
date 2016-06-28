function [ H, dH ] = Daniel_get_W_Hamiltonian( DATA, k )
% Erzeugt Hamiltonian und die Gradieten nach k_x, k_y. Benötigt einen
% k-Punkt, sowie die Struktur DATA, die die Felder W_DATA, die
% kartesische Basis und den degeneracy-Vektor enthalten muss.
%
% DATA:
%       DATA.W_DATA, DATA.RBASE, DATA.degeneracy

nOrbs = 3;

% Wannier-Daten
W_DATA     = DATA.W_DATA;
R_BASE     = DATA.RBASE;
degeneracy = DATA.degeneracy;

% Koordinatentransformation
TRF   = DATA.TRF;
k_TRF = TRF * k(:)
k_TRF = k_TRF(:).';

idx_deg = 1;

H =  zeros( nOrbs, nOrbs ) ;
dH =  zeros( 2*nOrbs, nOrbs ) ;

% Hamiltonian zusammensetzen
for ii = 1:1:size( W_DATA, 1 )
    
    % H(i,j) = H(i,j) + t(i,j) * exp( 2 pi i k*r(i,j) ) / deg(i,j)
    H_ij = W_DATA( ii, 6 ) * exp( 2 * pi * 1i * sum( k_TRF( 1:3 ) .* W_DATA( ii, 1:3 ) ) ) / degeneracy( idx_deg );
    
    H( W_DATA( ii, 4 ), W_DATA( ii, 5 ) ) = H( W_DATA( ii, 4 ), W_DATA( ii, 5 ) ) + H_ij;
    
    % d/dkx H = H(i,j) + i r(1) t(i,j) * exp( 2 pi i k*r(i,j) ) / deg(i,j)
    dH( W_DATA( ii, 4 ), W_DATA( ii, 5 ) ) = ...
        dH( W_DATA( ii, 4 ), W_DATA( ii, 5 ) ) + 1i * ( W_DATA( ii, 1 ) * R_BASE( 1, 1 ) + W_DATA( ii, 2 ) * R_BASE( 2, 1 ) ) * H_ij;
    
    % d/dky H = H(i,j) + i r(2) t(i,j) * exp( 2 pi i k*r(i,j) ) / deg(i,j)
    dH( W_DATA( ii, 4 ) + nOrbs, W_DATA( ii, 5 ) ) = ...
        dH( W_DATA( ii, 4 ) + nOrbs, W_DATA( ii, 5 ) ) + 1i * ( W_DATA( ii, 1 ) * R_BASE( 1, 2 ) + W_DATA( ii, 2 ) * R_BASE( 2, 2 ) ) * H_ij;
    
    tmp = ii / size( H, 1 )^2;
    
    if ( floor( tmp ) - ceil( tmp ) == 0 )
        
        idx_deg = idx_deg + 1;
        
    end
    
end

end