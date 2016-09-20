function [ Bloch , energy , nr_w ] = call_bloch_structure( Ctrl, constAg , Para , Data )

Bloch.hbar = constAg.hbar;                  % Planck

Bloch.E = Ctrl.E;                           % All information for E-field
Bloch.gamma = Ctrl.ode.dephrasing;          % Dephrasing

Bloch.nr.k = Para.nr.k;                     % # k points
Bloch.nr.d = Para.nr.dipol;                 % # dipole transitions
Bloch.ind = reshape( 1 : Para.nr.dipol * Para.nr.k ,[], Para.nr.dipol);     % Indexing of psi per transition

% Weight of k-points as coulumn vector
Bloch.wk = repmat( Data.wk.' * Para.BZsmall.area , [Para.nr.dipol,1] );     % in BZred for P
Bloch.wkentire = Data.wk.' * Para.BZsmall.area / 6;                         % in BZ for E-Rabi      - maybe out




Bloch.Esum = zeros(Para.nr.k * Para.nr.dipol , 1 );
% Bloch.feff = zeros(Para.nr.k * Para.nr.dipol , 1 );                         % In the linear regime. feff const.

for ii = 1:Para.nr.dipol
    Bloch.Esum( Bloch.ind(:,ii) ) = ( Data.Ek( Para.dipol_trans(ii,1),: , 1 ) + Data.Ek( Para.dipol_trans(ii,2),: , 1 ) ).' ;        
%     Bloch.feff( Bloch.ind(:,ii) ) = 1 - ( Data.fk(Para.dipol_trans(ii,1),:) + Data.fk(Para.dipol_trans(ii,2),:) ).';       
end

Dipol_array = cell2mat(Data.dipol);
Bloch.dipol = Dipol_array(:,Ctrl.E.polarization);

Bloch.coul_ctrl = Ctrl.Coul.active;
Bloch.coulomb = Data.V_rabi;
Bloch.coul_dip_mapping = Para.coul_dip_mapping;


energy = linspace( Ctrl.ode.energy_range(1) , Ctrl.ode.energy_range(2) , Ctrl.ode.energy_steps )';
Bloch.w = energy / constAg.hbar;                                            % Energiefenster in omega ???
Bloch.nr.w = numel(Bloch.w);
nr_w = Bloch.nr.w;
