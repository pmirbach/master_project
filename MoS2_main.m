%% Header
clear variables
clear global
% close all

% profile off
% profile on
% clc
dbstop if error

% Unterordner fuer Funktionsaufrufe:
addpath(genpath(pwd))

load('KonstantenAg.mat')    % Naturkonstanten (Ag Jahnke)

%% Steuerungsdatei

% Allgemeines
Ctrl.material = 'MoS2';         % Materialien: MoS2 WS2 MoSe2 WSe2 MoTe2 WTe2 
Ctrl.lattice_constant = 0.318;  % ab initio TB: MoS2: .316, .318, .320
% Calculations with Liu TB modell:  .319    !!


% k-mesh % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctrl.k_mesh_mp.qr = 30;        % Unterteilungsgröße
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
Ctrl.plot.tb = [1 , 1];         % Bandstructure
Ctrl.plot.exc = [0 , 0];         % Excitation
Ctrl.plot.dipol = [0 , 0];      % Dipol matrix elements
Ctrl.plot.coul = 0;
Ctrl.plot.ren_bs = [0 , 0];      % Dipol matrix elements

Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 0;         % 1 ganze BZ, 0 nur red. BZ

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

Para = call_para(Ctrl, constAg);

constAg.hbar = 0.6582119514;
% load('Pfad.mat')
% load('Pfad_35.mat')

%% Monkhorst-Pack

[ Data.k , Data.wk ] = k_mesh_mp(Ctrl, Para);
Para.nr.k = size(Data.k,2);

Para.symm_indices = find( Data.wk == 1 );

Data.k = round(Data.k,13);



%%



% kround = reshape( Data.k , 1, [] );
% for ii = 1:size(kround,2)
% %     disp(ii)
%     kstr = num2str( kround(ii) );
%     kstr(end) = '0';
%     kround(ii) = str2double( kstr );
% end
% knew = reshape( kround, size( Data.k ) );
% 
% Data.k = knew;


% % load kpts_11x11.mat
% % k1 = permute(k_11x11,[2,1,3]);
% % k1 = k1(:,1:66);
% 
% load kpts_35x35.mat
% k1 = permute(kpts_35x35,[2,1,3]);
% 
% k2 = k1(1:2,:);
% Data.wk = round( k1(3,:) / min(k1(3,:)) );
% Para.BZsmall.area = 1;
% Para.symm_indices = find( Data.wk == 1 );
% 
% [Data.k] = red_to_BZ(k2);
% Para.nr.k = size(Data.k,2);
% 
% % load kpts_11x11(2).mat
% % k3 = permute(kpts,[2,1,3]);
% % Data.k = k3(1:2,:,:);
% % Para.nr.k = size(Data.k,2);
% 
% Para.BZsmall.area = min(k1(3,:));
% % Para.BZsmall.area = min(k3(3,:));
% Para.coul.pol = 1.27287195103197  / Para.BZsmall.area;        % 35 x 35
% % Para.coul.pol = 4.32776463350871  / Para.BZsmall.area;          % 11 x 11
%  
% Para.k.qmin = 0.193102996717152;                                % 35 x 35
% % Para.k.qmin = 0.656550188837845;                                % 11 x 11


%%

kx = Data.k(1,:,1);
ky = Data.k(2,:,1);

%% Tight-Binding
fprintf('Tight-binding:        Start'); tic

[Data.Ek, Data.Ev, Prep.Ek_noSOC, Prep.Ev_noSOC] = tight_binding_liu(Ctrl, Para, Data);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Data.Ek(:,:,1),[2 3]);

%% Tight-Binding - ab initio

[Ek, Ev, Ek_noSOC, Ev_noSOC] = tight_binding_roesner(Ctrl, Para, Data);

%%
[fig.bandstr_surf2, fig.bandstr_path2] = plot_bandstr(Ctrl,Para,Data.k,Ek(:,:,1),[2 3]);


%%

% nEv = EV_ortho( Para, Data.k, Prep.Ek_noSOC , Prep.Ev_noSOC );

% Prep.Ev_noSOC = nEv;

%% Daniels Eigenvektoren
% load('CVec (2).mat')
% Data.Ev = CVec;
% Data.Ev = abs(real(Data.Ev)) + 1i * abs(imag(Data.Ev));

% load('CVec_35.mat')
% Data.Ev = CVec;

% load test_eig_chol.mat
% Data.Ek = ENERGY;
% Data.Ev = coeff;

% load('CVec_35_noSOC.mat')
% 
% compa = Ev_comp(Prep.Ev_noSOC, CVec(1:3,1:3,:,:));
% 
% figure
% set(gcf, 'Color', 'w');
% for ii = 1:3
%     subplot(1,3,ii)
%     imagesc(squeeze(sum(real(compa(ii,:,:)),3)))
%     colorbar
% end

% figure
% set(gcf, 'Color', 'w');
% for ii = 1:6
%     subplot(2,3,ii)
%     scatter3(kx,ky,imag(Data.Ev(1,2,:,ii)))
% end



% load('CVec_35_noSOC.mat')

% Prep.Ev_noSOC(:,1,:,:) = - Prep.Ev_noSOC(:,1,:,:);


%% Simulation-preperations
fprintf('Preperations:         Start'); tic

[Prep.Eks, Prep.CV, Prep.CV_noSOC, Prep.minq] = prep(Para, Data, Prep.Ev_noSOC);

fprintf('   -   Finished in %g seconds\n',toc)

%%
% Prep.CV = ones(size(Prep.CV));
% Prep.CV = call_CV( Ev )


%% Thermische Anregung
fprintf('Excitation:           Start'); tic

[Data.fk, Para.mu] = excitation(Ctrl,constAg,Para,Data.wk,Prep.Eks);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.exc_surf, fig.exc_path] = plot_excitation(Ctrl,Para,Data.k,Data.fk,[2 3]);


%% Coulomb WW
fprintf('Coulomb matrix:       Start'); tic

[Data.V.f, Data.V.h, Data.V.h_off] = coulomb_hf( Ctrl , Para , Prep );

fprintf('   -   Finished in %g seconds\n',toc)

%% Coulomb Plots

[fig.coulomb_up, fig.coulomb_down] = plot_coulomb( Ctrl , Para,Data.k , Data.V.f , Para.symm_indices(2) );

%% Band renorm

fprintf('Band renormalization: Start'); tic

[Data.Ek_hf, Data.Ek_h, Data.Ek_f, Test.ren_h] = renorm2(Para, Data.Ek, Data.V, Data.fk, Data.wk);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.ren_bandstr_surf, fig.ren_bandstr_path] = plot_renorm_bandstr(Ctrl,Para,Data.k,[Data.Ek(:,:,1);Data.Ek_h],[2 3]);

% as1 = plot_path(Ctrl,Para,Data.k,Test.ren_h,200);      % Test der Hartree Renormierung

%%
% as2 = plot_path(Ctrl,Para,Data.k,Prep.Eks,200);

%% structure für Variablen für Blochgleichungen

% Para.dipol_trans = [1, 2 ; 1 , 3 ; 4 , 5 ; 4 , 6 ];
Para.dipol_trans = [1 2; 4 5];
% Para.dipol_trans = [4 5];
Para.nr.dipol = size(Para.dipol_trans,1);


%% Dipolmatrix
fprintf('Dipol:                Start'); tic

Data.dipol = dipol(Para, Prep, Data);

fprintf('   -   Finished in %g seconds\n',toc)

titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol,[2 2],titlestr);
clear titlestr

%%

Bloch.nrd = Para.nr.dipol;
Bloch.ind = reshape( 1 : Para.nr.dipol*Para.nr.k ,[], Para.nr.dipol);


[V_rabi_fock] = coulomb_rabi_f(Ctrl, Para, Prep, Data.Ev);                           % All coulomb matrices. (in 3rd dimension)
Bloch.coulomb = V_rabi_fock;

% Hab ich schon
Bloch.hbar = constAg.hbar;
Bloch.wk = repmat( Data.wk.' * Para.BZsmall.area , [Para.nr.dipol,1] );     % Spaltenvektor
Bloch.wkentire = Data.wk.' * Para.BZsmall.area / 6;                         % Spaltenvektor
% Bloch.wkentire = Bloch.wk / 6; 



Bloch.Esum = zeros(Para.nr.k * Para.nr.dipol , 1 );
Bloch.dipol = zeros(Para.nr.k * Para.nr.dipol , 1 );
Bloch.feff = zeros(Para.nr.k * Para.nr.dipol , 1 );                         % In the linear regime. feff const.

for ii = 1:Para.nr.dipol
    Bloch.Esum( Bloch.ind(:,ii) ) = ( Prep.Eks( Para.dipol_trans(ii,1),: ) + Prep.Eks( Para.dipol_trans(ii,2),: ) ).';
    
%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{dipolnr}(1,:) - 1i * Data.dipol{dipolnr}(2,:) ).';
    Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).'; 
%     Bloch.dipol( (ii-1) * Para.nr.k + 1 : ii * Para.nr.k ) = abs( Data.dipol{ii}(1,:) ).'; 
    
    Bloch.feff( Bloch.ind(:,ii) ) = 1 - ( Data.fk(Para.dipol_trans(ii,1),:) + Data.fk(Para.dipol_trans(ii,2),:) ).';
    
end

% Bloch.dipol = 5e4 * ones(Para.nr.k * Para.nr.dipol,1);                % ? 1-3 too strong.


Bloch.gamma = 10;
Bloch.E0 = 1e-7;
Bloch.t_peak = 0.003;
Bloch.sigma = 0.001;
Bloch.nrk = Para.nr.k;



% Kommt noch dazu
Emin = -1000;
Emax = 0;
E = linspace(Emin,Emax,2001)';

Bloch.w = E / constAg.hbar;             % Energiefenster in omega ???

Para.nr.w = numel(Bloch.w);


%% Zeitentwicklung

Bloch.coul_ctrl = 1;                    % Coulomb Interaktion


tspan = [0 0.4];
% tspan = linspace(0, 0.4, 30000);
psik_E_ini = zeros(1 , Para.nr.dipol * Para.nr.k + size(Bloch.w,1) * 2);

% options=odeset('OutputFcn',@odeprog,'Events',@odeabort);
options=odeset('OutputFcn',@odeprog,'Events',@odeabort,'RelTol',1e-6,'AbsTol',1e-8);
% opts = odeset('RelTol',1e-1,'AbsTol',1e-3);
[t,psik_E] = ode113(@(t,psik_E) dgl_bloch(t,psik_E,Bloch), tspan, psik_E_ini, options);


%

% psik = psik_E(:,1:Parameter.nrk);
P_w = psik_E(end,( end - 2 * Para.nr.w + 1 ):( end - 1 * Para.nr.w ));
E_w = psik_E(end,( end - 1 * Para.nr.w + 1 ):end);

chi_w = P_w ./ E_w;

% close all

% figure
plot(E , imag(chi_w))

hold on
% load('spec_V_dip.mat')
% load ist_egal.mat
% plot(spec_12(:,1)*Bloch.hbar , spec_12(:,3),'r--')
% plot(spec(:,1)*Bloch.hbar , spec(:,3))


% xlim([-500 0])

%%
% P_t2 = zeros(1,numel(t));
% 
% for ii = 1:numel(t)
%     P_t2(ii) = 1 / (2 * pi)^2 * Bloch.wk.' * (conj(Bloch.dipol) .*  psik_E(ii,1:Bloch.nrd*Bloch.nrk).' );
% end
% 
% figure
% plot(real(P_t2))
% hold on
% plot(imag(P_t2))

% E_t2 = Bloch.E0 * exp(-1/2 * ( ( t.' - Bloch.t_peak ) / Bloch.sigma ).^2 * 4 * log(2) );
%
% P_w2 = fft(P_t2);
% E_w2 = fft(E_t2);
%
% chi_w2 = P_w2 ./ E_w2;
%
% close all
% plot(imag(P_w2))

%%

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
