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

load('KonstantenAg.mat')    % Naturkonstanten (Ag Jahnke)

%% Steuerungsdatei

% Allgemeines
Ctrl.material = 'MoS2';   % Materialien: MoS2 WS2 MoSe2 WSe2 MoTe2 WTe2
Ctrl.method = 'TNN';      % Möglich:   NN , TNN
Ctrl.SOC = 1;             % Spin-Orbit-Coupling

% k-mesh % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctrl.k_mesh_mp.qr = 18;        % Unterteilungsgröße
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% muss durch 6 teilbar sein, damit Hochsymmetriepunkte mit im mesh sind
% 60 -> 631 kpts; 120 -> 2461 kpts

% Anregungsdichte
Ctrl.temperature = 300;         % Temperatur in K
% Ctrl.carrier_density = 1e13;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density = 0;    % Anregungsdichte in 1/cm^2
Ctrl.carrier_density_tol = Ctrl.carrier_density * 1e-8;

%% Plot Control

Ctrl.plot.path = {'K','M', 'K*', '\Gamma', 'K','M','\Gamma'};
% Ctrl.plot.path = {'K' '\Gamma' 'K*'};

Ctrl.plot.k_mesh = [0 , 0];     % Kontrollbilder
% 1: Surface, 2: Pathplot
Ctrl.plot.tb = [0 , 1];         % Bandstructure
Ctrl.plot.exc = [0 , 0];         % Excitation
Ctrl.plot.dipol = [0 , 0];      % Dipol matrix elements
Ctrl.plot.coul = 0;
Ctrl.plot.ren_bs = [0 , 0];      % Dipol matrix elements

Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 0;         % 1 ganze BZ, 0 nur red. BZ

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

Para = call_para(Ctrl, constAg);

%% Monkhorst-Pack

% [Data.k, Data.wk] = k_mesh_mp(Ctrl, Para);
% Para.nr.k = size(Data.k,2);
% 
% Para.symm_indices = find( Data.wk == 1 );


load kpts_11x11.mat
k1 = permute(k_11x11,[2,1,3]);
k2 = k1(1:2,1:66);
Data.wk = k1(3,1:66);
Para.BZsmall.area = 1;
Para.symm_indices = find( Data.wk == 1 );

[Data.k] = red_to_BZ(k2);
Para.nr.k = size(Data.k,2);

 
% Para.BZsmall.area = (8 / 3 / sqrt(3) / ( 55 - 1 )^2)*(pi/0.319)^2;
% Para.coul.pol = 0.801437895090000 / Para.BZsmall.area;
% Para.k.qmin = 0.1216;


%% Tight-Binding
fprintf('Tight-binding:        Start'); tic

[Data.Ek, Data.Ev, Prep.Ek_noSOC, Prep.Ev_noSOC] = tight_binding_liu(Ctrl, Para, Data);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Data.Ek(:,:,1),[2 3]);

%% Simulation-preperations
fprintf('Preperations:         Start'); tic

[Prep.Eks, Prep.CV, Prep.CV_noSOC, Prep.minq] = prep(Para, Data, Prep.Ev_noSOC);

fprintf('   -   Finished in %g seconds\n',toc)

%% Thermische Anregung
fprintf('Excitation:           Start'); tic

[Data.fk, Para.mu] = excitation(Ctrl,constAg,Para,Data.wk,Prep.Eks);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.exc_surf, fig.exc_path] = plot_excitation(Ctrl,Para,Data.k,Data.fk,[2 3]);

%% Dipolmatrix
fprintf('Dipol:                Start'); tic

Data.dipol = dipol(Para, Prep, Data);

fprintf('   -   Finished in %g seconds\n',toc)

titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','2 \rightarrow 1 \uparrow','2 \rightarrow 1 \uparrow'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol(1:3,1:3),[2 2],titlestr);
titlestr = {'1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow','2 \rightarrow 1 \downarrow','2 \rightarrow 1 \downarrow'};
[fig.dipolDown_surf, fig.dipolDown_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol(4:6,4:6),[2 2],titlestr);
clear titlestr

%% Coulomb WW
fprintf('Coulomb matrix:       Start'); tic

[Data.V.f, Data.V.h, Data.V.h_off] = coulomb_hf( Ctrl , Para , Prep );

fprintf('   -   Finished in %g seconds\n',toc)


%% Coulomb Plots

% [fig.coulomb_up, fig.coulomb_down] = plot_coulomb( Ctrl , Para,Data.k , Data.V.h , Para.symm_indices(2) );


%% Band renorm

fprintf('Band renormalization: Start'); tic

[Data.Ek_hf, Data.Ek_h, Data.Ek_f, Test.ren_h] = renorm2(Para, Data.Ek, Data.V, Data.fk, Data.wk);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.ren_bandstr_surf, fig.ren_bandstr_path] = plot_renorm_bandstr(Ctrl,Para,Data.k,[Data.Ek(:,:,1);Data.Ek_h],[2 3]);

% as1 = plot_path(Ctrl,Para,Data.k,Test.ren_h,200);      % Test der Hartree Renormierung


%%
as2 = plot_path(Ctrl,Para,Data.k,Prep.Eks,200);

%% structure für Variablen für Blochgleichungen
% Hab ich schon
Bloch.hbar = constAg.hbar;
Bloch.wk = Data.wk;

Bloch.Eks = Prep.Eks;

Bloch.dipol = 1 / sqrt(2) * abs( Data.dipol{2,1}(1,:) - 1i * Data.dipol{2,1}(2,:) ).';
% Bloch.dipol = 5e4 * ones(2461,1);


Bloch.gamma = 10;
Bloch.E0 = 1e-6;
Bloch.t_peak = 0.003;
Bloch.sigma = 0.001;
Bloch.nrk = Para.nr.k;



% Kommt noch dazu
Emin = -200;
Emax = 1000;
E = linspace(Emin,Emax,500)';

Bloch.w = E / constAg.hbar;             % Energiefenster in omega ???

Para.nr.w = numel(Bloch.w);


% Zeitentwicklung

tspan = [0 2];
psik_E_ini = zeros(1,Para.nr.k + size(Bloch.w,1) * 2);

options=odeset('OutputFcn',@odeprog,'Events',@odeabort);
% opts = odeset('RelTol',1e-1,'AbsTol',1e-3);
[t,psik_E] = ode45(@(t,psik_E) dgl_bloch(t,psik_E,Bloch), tspan, psik_E_ini, options);


%

% psik = psik_E(:,1:Parameter.nrk);
P_w = psik_E(end,( end - 2 * Para.nr.w + 1 ):( end - 1 * Para.nr.w ));
E_w = psik_E(end,( end - 1 * Para.nr.w + 1 ):end);

chi_w = P_w ./ E_w;

close all
plot(E * 1e-3, imag(chi_w))

%%
P_t2 = zeros(1,numel(t));
% psik = psik_E(1:Bloch.nrk);

for ii = 1:numel(t)
    P_t2(ii) = 1 / (2 * pi)^2 * Data.wk * (conj(Bloch.dipol) .*  psik_E(ii,1:Bloch.nrk).' );
end

E_t2 = Bloch.E0 * exp(-1/2 * ( ( t.' - Bloch.t_peak ) / Bloch.sigma ).^2 * 4 * log(2) );

P_w2 = fft(P_t2);
E_w2 = fft(E_t2);

chi_w2 = P_w2 ./ E_w2;

close all
plot(imag(P_w2))

% figure
% plot(t,real(P_t2))
% hold on
% plot(t,imag(P_t2),'r')
% legend('real','imag')

%%
% plot_surf(Para,Data.k(:,:,1),Bloch.dipol,[1,1])


%%
% profile viewer
% profile off
