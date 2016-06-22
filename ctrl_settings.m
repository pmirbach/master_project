function Ctrl = ctrl_settings


%%%%%%%%%%%%% Material - lattice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
Ctrl.material = 'MoS2';         % Materialien: MoS2 WS2 MoSe2 WSe2 MoTe2 WTe2 
Ctrl.lattice_constant = 0.318;  % ab initio TB: MoS2: .316, .318, .320   (in liu a = 0.319 auto)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% k-mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
Ctrl.k_mesh_mp.qr = 60;        % Unterteilungsgröße (Vielfaches von 6!). 60 -> 631 kpts; 120 -> 2461 kpts

Ctrl.comp_Daniel.k = 0;
Ctrl.comp_daniel.Ev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Tight-Binding modell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
Ctrl.TB_modell = 'ab_initio';   % Tight-Binding modell: Roesner ab initio or liu
% Ctrl.TB_modell = 'liu';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Coulomb interaction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
Ctrl.eps_r = 1;                 % Permiabilitaet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anregungsdichte
Ctrl.temperature = 300;         % Temperatur in K
% Ctrl.carrier_density = 1e13;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density = 0;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density_tol = Ctrl.carrier_density * 1e-8;



Ctrl.profile_flag = 0;

%% Plot Control

Ctrl.plot.path = {'K','M', 'K*', '\Gamma', 'K','M','\Gamma'};
% Ctrl.plot.path = {'K' '\Gamma' 'K*'};

Ctrl.plot.k_mesh = [0 , 0];     % Kontrollbilder
% 1: Surface, 2: Pathplot
Ctrl.plot.tb = [0 , 0];         % Bandstructure
Ctrl.plot.exc = [0 , 0];         % Excitation
Ctrl.plot.dipol = [1 , 0];      % Dipol matrix elements
Ctrl.plot.coul = 0;
Ctrl.plot.ren_bs = [0 , 0];      % Dipol matrix elements

Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 1;         % 1 ganze BZ, 0 nur red. BZ