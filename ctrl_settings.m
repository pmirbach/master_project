function Ctrl = ctrl_settings


%%%%%%%%%%%%% Material - lattice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% Ctrl.material = 'MoS2';         % Materialien: MoS2 WS2 MoSe2 WSe2 ( MoTe2 WTe2 later )
% Ctrl.material_trans_me = 'Mo';
% % ab initio TB: MoS2:   0.3160, 0.3180, 0.3200
% %               MoSe2:  0.3320
% %               WS2:    0.3191
% %               WSe2:   0.3325
% Ctrl.lattice_constant = 0.3180;

% Ctrl.material = 'MoS2';
% Ctrl.material_trans_me = 'Mo';
% Ctrl.lattice_constant = 0.3180;
% Ctrl.material_lambda = 0.074; 

Ctrl.material = 'MoSe2';
Ctrl.material_trans_me = 'Mo';
Ctrl.lattice_constant = 0.3320;
Ctrl.material_lambda = 0.232; 

% Ctrl.material = 'WS2';
% Ctrl.material_trans_me = 'W';
% Ctrl.lattice_constant = 0.3191;
% Ctrl.material_lambda = 0.2135; 

% Ctrl.material = 'WSe2';
% Ctrl.material_trans_me = 'W';
% Ctrl.lattice_constant = 0.3325;
% Ctrl.material_lambda = 0.232; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% k-mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Ctrl.k_mesh.type = 'symm';     % Orientation of k-mesh:   liu, symm
% Ctrl.k_mesh.type = 'liu';
Ctrl.k_mesh.qr = 90;        % Unterteilungsgröße (Vielfaches von 6!). 60 -> 631 kpts; 120 -> 2461 kpts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Tight-Binding modell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
Ctrl.TB_modell = 'ab_initio';   % Tight-Binding modell: Roesner ab initio or liu
% Ctrl.TB_modell = 'liu';
Ctrl.TB_t_symm = 1;             % Time symmetrization (only possible with k-mesh: symm) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Coulomb interaction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Ctrl.Coul.active = 1;                % Coulomb interaction active or not
Ctrl.Coul.eps_r = 1;                 % Permiabilitaet

Ctrl.Coul.Resta_fit_eps = 1;
Ctrl.Coul.eps_2 = 1;
Ctrl.Coul.eps_3 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anregungsdichte
Ctrl.temperature = 300;         % Temperatur in K
% Ctrl.carrier_density = 1e13;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density = 0;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density_tol = Ctrl.carrier_density * 1e-8;


Ctrl.cmp.load = 0;
Ctrl.cmp.use_k = 0;
% Ctrl.cmp.use_ev = 1;

Ctrl.profile_flag = 0;

%% Plot Control

Ctrl.plot.path = {'K','M', 'K*', '\Gamma', 'K','M','\Gamma'};
% Ctrl.plot.path = {'K' '\Gamma' 'K*'};

Ctrl.plot.k_mesh = [0 , 0];     % Kontrollbilder
% 1: Surface, 2: Pathplot
Ctrl.plot.tb = [0 , 0];         % Bandstructure
Ctrl.plot.exc = [0 , 0];         % Excitation
Ctrl.plot.dipol = [0 , 0];      % Dipol matrix elements
Ctrl.plot.coul = 0;
Ctrl.plot.ren_bs = [0 , 0];      % Dipol matrix elements

Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 0;         % 1 ganze BZ, 0 nur red. BZ