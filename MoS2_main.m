%% Header
clear variables
clear global
close all
% profile off
% profile on
% clc
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
Ctrl.k_mesh_mp.qr = 18;        % Unterteilungsgröße

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
Ctrl.plot.exc = [0 , 0];         % Excitation
Ctrl.plot.dipol = [0 , 0];      % Dipol matrix elements


Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 0;         % 1 ganze BZ, 0 nur red. BZ

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

[Parameter.TB.liu.names, Parameter.TB.liu.values] = load_TB_parameter(Ctrl);

Parameter.symmpts{1} = {'\Gamma', 'K', 'K*', 'M'};
Parameter.symmpts{2} = 2 * pi / (3 * Parameter.TB.liu.values(1)) ...
    * [0, 0 ; 2, 0 ; 1, sqrt(3) ; 3 / 2, sqrt(3) / 2 ]';
Parameter.rezGV = 2 * pi / Parameter.TB.liu.values(1) * [1, -1 / sqrt(3); 0, 2 / sqrt(3)]';
Parameter.area_real = 3 * sqrt(3) / 2 * (Parameter.TB.liu.values(1))^2;
Parameter.area_BZ = 3 * sqrt(3) / 2 * (norm(Parameter.symmpts{2}(:,2)))^2;
Parameter.area_sBZ = 3 * sqrt(3) / 2 * (norm(Parameter.symmpts{2}(:,2)) / Ctrl.k_mesh_mp.qr)^2;
Parameter.coul_screened = [[1.17 , 7.16 , 0.199, 2.675];
    [0.456 , 14.03 , 0.242, 2.930]; [0.46 , 13.97 , 0.24, 2.930]; 
    [1.288 , 6.88 , 0.232, 2.682]; [0.713 , 9.9 , 0.246, 2.855]; 
    [1.273 , 6.92 , 0.227, 2.682]];
Parameter.coul_kappa = 1000;             % Kappa, because of Singularity
Parameter.dipol_trans = [1, 2 ; 1 , 3 ; 2, 1 ; 3 , 1 ; 4 , 5 ; 4 , 6 ; 5 , 4 ; 6 , 4 ];

%% Monkhorst-Pack
[Data.k] = k_mesh_mp(Ctrl, Parameter);
Parameter.nrk = size(Data.k,2);

%% Tight-Binding
fprintf('Tight-binding:  Start'); tic

[Data.Ek, Data.Ev] = tight_binding_liu(Ctrl, Parameter, Data);

period = toc; fprintf('   -   Finished in %g seconds\n',period) 

[fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Parameter,Data.k,Data.Ek(:,:,1),[2 3]);

%% Simulation-preperations
fprintf('Preperations:   Start'); tic

[Prep.Eks, Prep.CV, Prep.minq, Prep.coul_intrp, Prep.coul_h] = prep(constAg, Parameter, Data);

period = toc; fprintf('   -   Finished in %g seconds\n',period) 


%% Thermische Anregung
fprintf('Excitation:     Start'); tic

Data.fk = excitation(Ctrl,constAg,Data.k(:,:,1),Prep.Eks);

period = toc; fprintf('   -   Finished in %g seconds\n',period) 

[fig.exc_surf, fig.exc_path] = plot_excitation(Ctrl,Parameter,Data.k,Data.fk,[2 3]);


%% Dipolmatrix
fprintf('Dipol:          Start'); tic

Data.dipol = dipol(Parameter, Prep, Data);

period = toc; fprintf('   -   Finished in %g seconds\n',period) 

titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','2 \rightarrow 1 \uparrow','2 \rightarrow 1 \uparrow'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Parameter,Data.k,Data.dipol(1:3,1:3),[2 2],titlestr);
titlestr = {'1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow','2 \rightarrow 1 \downarrow','2 \rightarrow 1 \downarrow'};
[fig.dipolDown_surf, fig.dipolDown_path] = plot_dipol(Ctrl,Parameter,Data.k,Data.dipol(4:6,4:6),[2 2],titlestr);


%% Tests
% figure; hold on
% jj = 2;
% for ii = 1:6
%     plot(Prep.minq(jj,:,ii))
% end

%% Coulomb WW

% profile on

% [Ek_hf,Ek_h,Ek_f] = coulomb_1(constAg,Parameter,Data,Prep.CV);
[Ek_hf,Ek_h,Ek_f] = coulomb_2(constAg,Parameter,Data,Prep);

% profile report
% profile off

%% 


% [d1, d2] = plot_bandstr(Ctrl,Parameter,Data.k,Prep.Eks,[2 3]);


% y0 = 0;

%% Zeitentwicklung

% tspan = [0 2];
% psik_E_ini = zeros(1,Parameter.nrk + 4002);
% 
% % opts = odeset('RelTol',1e-1,'AbsTol',1e-3);
% [t,psik_E] = ode45(@(t,psik_E) dgl_bloch(t,psik_E,Data,constAg), tspan, psik_E_ini);
% % plot(t,y,'-o')
 
%%

% psik = psik_E(:,1:Parameter.nrk);
% P_w = psik_E(end,Parameter.nrk+1:Parameter.nrk+2001);
% E_w = psik_E(end,Parameter.nrk+2002:end);
% 
% chi_w = P_w ./ E_w;
% 
% plot(-1000:1000,imag(chi_w))


% P = zeros(1,numel(t));
% 
% d = 1 / sqrt(2) * transpose(Data.dipol{2,1}(1,:) - 1i * Data.dipol{2,1}(2,:));
% 
% for ii = 1:numel(t)
%     P(ii) = 1 / (2 * pi)^2 * Data.k(3,:,1) * (conj(d) .* transpose(psik));
% end
% 
% figure
% plot(t,real(P))
% hold on
% plot(t,imag(P),'r')
% legend('real','imag')




%%
% profile viewer
% profile off
