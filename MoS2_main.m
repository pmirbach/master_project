%% Header
clear variables
clear global
close all
% profile on
clc
dbstop if error

% Unterordner für Funktionsaufrufe:
addpath(genpath(pwd))
% addpath(genpath('/home/pmirbach/Masterarbeit/M_bibliothek'))

load('KonstantenAg.mat')    % Naturkonstanten (Alex und Co)

%% Steuerungsdatei

% Allgemeines
Ctrl.material = 'MoS2';   % Materialien: MoS2 WS2 MoSe2 WSe2 MoTe2 WTe2
Ctrl.method = 'TNN';      % Möglich:   NN , TNN
Ctrl.SOC = 1;             % Spin-Orbit-Coupling

% k-mesh
Ctrl.k_mesh_mp.qr = 60;        % Unterteilungsgröße
% muss durch 6 teilbar sein, damit Hochsymmetriepunkte mit im mesh sind
% 60 -> 631 kpts; 120 -> 2461 kpts

% Anregungsdichte
Ctrl.temperature = 300;         % Temperatur in K
Ctrl.carrier_density = 1e13;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density_tol = Ctrl.carrier_density * 1e-8;

%% Plot Control

Ctrl.plot.path = {'\Gamma' 'K' 'M' 'K*' '\Gamma' 'M'};
% Ctrl.plot.path = {'K' '\Gamma' 'K*'};

Ctrl.plot.k_mesh = [0 , 0];     % Kontrollbilder
% 1: Surface, 2: Pathplot
Ctrl.plot.tb = [0 , 0];         % Bandstructure
Ctrl.plot.exc = [0, 0];         % Excitation
Ctrl.plot.dipol = [0 , 0];      % Dipol matrix elements


Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 0;         % 1 ganze BZ, 0 nur red. BZ

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

[Parameter.TB.liu.names, Parameter.TB.liu.values] = load_TB_parameter(Ctrl);

Parameter.symmpts{1} = {'\Gamma', 'K', 'K*', 'M'};
Parameter.symmpts{2} = 2 * pi / (3 * Parameter.TB.liu.values(1)) ...
    * [0, 0 ; 2, 0 ; 1, sqrt(3) ; 3 / 2, sqrt(3) / 2 ]';
Parameter.rezGV = 2 * pi / Parameter.TB.liu.values(1) ...
    * [1, -1 / sqrt(3); 0, 2 / sqrt(3)]';
Parameter.area_real = 3 * sqrt(3) / 2 * (Parameter.TB.liu.values(1))^2;
Parameter.area_BZ = 3 * sqrt(3) / 2 * (norm(Parameter.symmpts{2}(:,2)))^2;
Parameter.area_sBZ = 3 * sqrt(3) / 2 * ...
    (norm(Parameter.symmpts{2}(:,2)) / Ctrl.k_mesh_mp.qr)^2;
Parameter.coul_screened = [[1.17 , 7.16 , 0.199, 2.675];
    [0.456 , 14.03 , 0.242, 2.930]; [0.46 , 13.97 , 0.24, 2.930]; 
    [1.288 , 6.88 , 0.232, 2.682]; [0.713 , 9.9 , 0.246, 2.855]; 
    [1.273 , 6.92 , 0.227, 2.682]];

%% Monkhorst-Pack
[Data.k] = k_mesh_mp(Ctrl, Parameter);

%% Tight-Binding
[Data.Ek,Data.Ev] = tight_binding_liu(Ctrl, Parameter, Data);
[fig.bandstr_surf, fig.bandstr_path] = ...
    plot_bandstr(Ctrl,Parameter,Data.k,Data.Ek,[2 3]);

%% Thermische Anregung
Data.fk = excitation(Ctrl,constAg,Data.k(:,:,1),Data.Ek);
[fig.exc_surf, fig.exc_path] = ...
    plot_excitation(Ctrl,Parameter,Data.k,Data.fk,[2 3]);

%% Dipolmatrix

% Data.dipol = dipol(Parameter, Data);
% Darstellung
%% Darstellung Dipol
% A = 1 / sqrt(2) * (Data.dipol(1,:) + 1i * Data.dipol(2,:));
% BZ_surf(Ctrl,Data.k,abs(A))

% von = 1;
% zu = 2;
%
% A = 1 / sqrt(2) * (Data.dipol{von,zu}(1,:) + 1i * Data.dipol{von,zu}(2,:));
% BZ_surf(Ctrl,Data.k,abs(A),'normal')


% von = 2;
% zu = 3;
%
% A = 1 / sqrt(2) * (Data.dipol{von,zu}(1,:) + 1i * Data.dipol{von,zu}(2,:));
% BZ_surf(Ctrl,Data.k,abs(A))

% BZ_surf(Ctrl,Data.k,Data.dipol_plot(1:4,:),'subplot')
% BZ_surf(Ctrl,Data.k,Data.dipol_plot(5:8,:),'subplot')

% test= squeeze(P(3,3,:));

% BZ_surf(Ctrl,Data.k,abs(test),'at')

% profile report

%% Coulomb WW
coulomb_1(constAg,Parameter,Data)

%% Flächeninhalt
[B, B_integ] = flaecheninhalt(Parameter,Data.k(:,:,1));

%%
profile off
