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
Ctrl.k_mesh_mp.qr = 60;        % Unterteilungsgröße
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Ctrl.plot.dipol = [1 , 0];      % Dipol matrix elements
Ctrl.plot.ren_bs = [0 , 0];      % Dipol matrix elements

Ctrl.plot.save = 0;             % 1 Speichern, 0 nicht
Ctrl.plot.entireBZ = 0;         % 1 ganze BZ, 0 nur red. BZ

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

Para = call_para(Ctrl, constAg);

%% Monkhorst-Pack

[Data.k] = k_mesh_mp(Ctrl, Para);
Para.nr.k = size(Data.k,2);

%% Tight-Binding
fprintf('Tight-binding:  Start'); tic

[Data.Ek, Data.Ev, Prep.Ek_noSOC, Prep.Ev_noSOC] = tight_binding_liu(Ctrl, Para, Data);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Data.Ek(:,:,1),[2 3]);

%% Simulation-preperations
fprintf('Preperations:   Start'); tic

[Prep.Eks, Prep.CV, Prep.CV_noSOC, Prep.minq] = prep(Para, Data, Prep.Ev_noSOC);

fprintf('   -   Finished in %g seconds\n',toc)

%% Thermische Anregung
fprintf('Excitation:     Start'); tic

Data.fk = excitation(Ctrl,constAg,Para,Data.k(:,:,1),Prep.Eks);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.exc_surf, fig.exc_path] = plot_excitation(Ctrl,Para,Data.k,Data.fk,[2 3]);

%% Dipolmatrix
fprintf('Dipol:          Start'); tic

Data.dipol = dipol(Para, Prep, Data);

fprintf('   -   Finished in %g seconds\n',toc)

titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','2 \rightarrow 1 \uparrow','2 \rightarrow 1 \uparrow'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol(1:3,1:3),[2 2],titlestr);
titlestr = {'1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow','2 \rightarrow 1 \downarrow','2 \rightarrow 1 \downarrow'};
[fig.dipolDown_surf, fig.dipolDown_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol(4:6,4:6),[2 2],titlestr);

%% Coulomb WW
fprintf('Coulomb matrix: Start'); tic

[V_fock, V_hartree] = coulomb_5(constAg,Para,Data,Prep);

fprintf('   -   Finished in %g seconds\n',toc)


%% Band renorm
[A,B] = meshgrid(1:3,1:3);
c=cat(2,A,B);
ll=reshape(c,[],2);
ll = [ll; ll+3];

ll = ll(:,:);

ll2 = [ll, ll(:,2)];
ll2(ll2(:,3)>3,3) = ll2(ll2(:,3)>3,3)-3;

[Ek_hf, Ek_h, Ek_f] = renorm2(Para, Data.Ek, V_fock, V_hartree, Data.fk, Data.k(3,:,1),ll2);

close all
[fig.ren_bandstr_surf, fig.ren_bandstr_path] = plot_renorm_bandstr(Ctrl,Para,Data.k,[Data.Ek(:,:,1);Ek_f],[2 3]);



%%
[A,B] = meshgrid(1:3,1:3);
c=cat(2,A,B);
ll=reshape(c,[],2);
ll = [ll; ll+3];

d = 0;

figure
for ii = 1:9
    subplot(3,3,ii)
    scatter3(Data.k(1,:,1),Data.k(2,:,1),V_fock(1,:,ii)')
    hold on
    scatter3(Data.k(1,:,1),Data.k(2,:,1),V_fock(1,:,ii+9)','+')
    title(num2str(ll(ii+d,:)))
end




%%
[B,B_integ] = flaecheninhalt(Para,Data.k);

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
