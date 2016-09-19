function [ Bloch ] = call_bloch_structure( Ctrl, constAg , Para , Data )

Bloch.nrd = Para.nr.dipol;
Bloch.ind = reshape( 1 : Para.nr.dipol * Para.nr.k ,[], Para.nr.dipol);

Bloch.coulomb = V_rabi;

% Hab ich schon
Bloch.hbar = constAg.hbar;
Bloch.wk = repmat( Data.wk.' * Para.BZsmall.area , [Para.nr.dipol,1] );     % Spaltenvektor
Bloch.wkentire = Data.wk.' * Para.BZsmall.area / 6;                         % Spaltenvektor
% Bloch.wkentire = Bloch.wk / 6; 



Bloch.Esum = zeros(Para.nr.k * Para.nr.dipol , 1 );
% Bloch.dipol = zeros(Para.nr.k * Para.nr.dipol , 1 );
Bloch.feff = zeros(Para.nr.k * Para.nr.dipol , 1 );                         % In the linear regime. feff const.


Bloch.dipol = A(:,1);

Bloch.coul_dip_mapping = Para.coul_dip_mapping;

for ii = 1:Para.nr.dipol
    Bloch.Esum( Bloch.ind(:,ii) ) = ( Data.Ek( Para.dipol_trans(ii,1),: , 1 ) + Data.Ek( Para.dipol_trans(ii,2),: , 1 ) ).' ;        
%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).'; 

%     Rausgek�rzt V
%     Bloch.feff( Bloch.ind(:,ii) ) = 1 - ( Data.fk(Para.dipol_trans(ii,1),:) + Data.fk(Para.dipol_trans(ii,2),:) ).';      
    
end


Bloch.gamma = 10;           % Dephrasing


Bloch.E0 = 1e-7;
Bloch.t_peak = 2.038 * 1e-3;
Bloch.sigma = 1e-3;
Bloch.nrk = Para.nr.k;