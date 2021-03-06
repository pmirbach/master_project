%% Header
clear variables
clear global
% close all


profile off
dbstop if error

% Unterordner fuer Funktionsaufrufe:
addpath(genpath(pwd))


%% Steuerungsdatei

Ctrl = ctrl_settings;

if Ctrl.profile_flag == 1
    profile on
end

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

% Einheitensystem:  mEv,  Ang,  ps,  pA


load('KonstantenAg.mat')    % Naturkonstanten (Ag Jahnke)
[ Para , W90Data , Coul_ME ] = call_para(Ctrl, constAg);

% Achtung: orbital order in Maltes TB modell different
% Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);                                             % ??? Kein Unterschied???






%% k - mesh
q_max = 2;
N_q = 100;
N_phi = 30;

Ctrl.k_mesh.qr = 240;

[ Data.k , Data.wk , Para.nr.k , index_isoc ] = k_mesh_isotrop( q_max, N_q, N_phi );

% plot(Data.k(1,:),Data.k(2,:),'x')



%%


% file = 'Daniel/MoS2_gradtest.mat';
% [Para, Data, Daniel] = load_daniel_data(Ctrl, Para, Data, file); % overwrites k, wk, coul_pol, wk_area, nr_k
% 
% % Fast kx, ky for scatter3 plots:
kx = Data.k(1,:);
ky = Data.k(2,:);




%% Tight-Binding
Ctrl.TB_t_symm = 0;

fprintf('Tight-binding (ab initio): Start'); tic

[Data.Ek, Data.Ev, Data.EGap, Prep.Ev_noSOC, Prep.H_grad_kx, Prep.H_grad_ky, Ek_old] = tight_binding_roesner_1d(Ctrl, Para, Data.k, W90Data);

fprintf('   -   Finished in %g seconds\n',toc)

% [fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Data.Ek(:,:,1),[2 3]);
% [fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Ek_old(:,:,1),[2 3]);

%%
% scatter3( Data.k(1,:), Data.k(2,:) , Data.Ek(1,:) )


%% Simulation-preperations
fprintf('Preperations:              Start'); tic

[Prep.CV, Prep.CV_noSOC, Prep.minq, Prep.V_ab_interpl ] = prep( Ctrl, Para, Data, Prep.Ev_noSOC , Coul_ME );

fprintf('   -   Finished in %g seconds\n',toc)

% [fig.bandstr_surf2, fig.bandstr_path2] = plot_bandstr(Ctrl,Para,Data.k,Prep.Eks,[1 2]);


%% Dipolmatrix
fprintf('Dipol:                     Start'); tic

Data.dipol = dipol(Para, Prep, Data);

fprintf('   -   Finished in %g seconds\n',toc)

titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol,[2 2],titlestr);
clear titlestr

Ploter.dipol = zeros( Para.nr.k, Para.nr.dipol );
for ii = 1:Para.nr.dipol
    Ploter.dipol_l(:,ii) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) + 1i * Data.dipol{ii}(2,:) ).'; 
    Ploter.dipol_r(:,ii) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).'; 
end

%%
1

%% Anregung, Bandrenormierung (Coulomb)     -   keine Zeit mehr
% %% Thermische Anregung
% fprintf('Excitation:                Start'); tic
% 
% [Data.fk, Para.mu] = excitation(Ctrl,constAg,Para,Data.wk,Prep.Eks);
% 
% fprintf('   -   Finished in %g seconds\n',toc)
% 
% [fig.exc_surf, fig.exc_path] = plot_excitation(Ctrl,Para,Data.k,Data.fk,[2 3]);
% 
% 
% %% Coulomb WW
% fprintf('Coulomb matrix:            Start'); tic
% 
% [Data.V.f, Data.V.h, Data.V.h_off] = coulomb_hf( Ctrl , Para , Prep );
% 
% fprintf('   -   Finished in %g seconds\n',toc)
% 
% %% Coulomb Plots
% 
% [fig.coulomb_up, fig.coulomb_down] = plot_coulomb( Ctrl , Para,Data.k , Data.V.f , Para.k_ind.symm(2) );
% 
% %% Band renorm
% 
% fprintf('Band renormalization:      Start'); tic
% 
% [Data.Ek_hf, Data.Ek_h, Data.Ek_f, Test.ren_h] = renorm2(Para, Data.Ek, Data.V, Data.fk, Data.wk);
% 
% fprintf('   -   Finished in %g seconds\n',toc)
% 
% [fig.ren_bandstr_surf, fig.ren_bandstr_path] = plot_renorm_bandstr(Ctrl,Para,Data.k,[Data.Ek(:,:,1);Data.Ek_h],[2 3]);
% 
% % as1 = plot_path(Ctrl,Para,Data.k,Test.ren_h,200);      % Test der Hartree Renormierung


%% structure für Variablen für Blochgleichungen

Bloch.nrd = Para.nr.dipol;
Bloch.ind = reshape( 1 : Para.nr.dipol*Para.nr.k ,[], Para.nr.dipol);

% tic
% [V_rabi_fock] = coulomb_rabi_f(Ctrl, Para, Prep, Data.Ev);                           % All coulomb matrices. (in 3rd dimension)
% toc
tic
[V_rabi_fock_interpl] = coulomb_rabi_f_interpl(Ctrl, Para, Prep );
toc
% %%
% figure; hold on
% scatter3(kx,ky,V_rabi_fock_interpl(1,:,1) )
% scatter3(kx,ky,V_rabi_fock(1,:,1),'r')
% % % figure; scatter3(kx,ky,V_rabi_fock_interpl(1,:,1) )

%%

Bloch.coulomb = V_rabi_fock_interpl;
% Bloch.coulomb = V_rabi_fock;

% Hab ich schon
Bloch.hbar = constAg.hbar;
Bloch.wk = repmat( Data.wk.' * Para.BZsmall.area , [Para.nr.dipol,1] );     % Spaltenvektor
Bloch.wkentire = Data.wk.' * Para.BZsmall.area / 6;                         % Spaltenvektor
% Bloch.wkentire = Bloch.wk / 6; 



Bloch.Esum = zeros(Para.nr.k * Para.nr.dipol , 1 );
Bloch.dipol = zeros(Para.nr.k * Para.nr.dipol , 1 );
Bloch.feff = zeros(Para.nr.k * Para.nr.dipol , 1 );                         % In the linear regime. feff const.

for ii = 1:Para.nr.dipol
    Bloch.Esum( Bloch.ind(:,ii) ) = ( Data.Ek( Para.dipol_trans(ii,1),: , 1 ) + Data.Ek( Para.dipol_trans(ii,2),: , 1 ) ).' ;
%     Bloch.Esum( Bloch.ind(:,ii) ) = ( - Data.Ek( Para.dipol_trans(ii,1),: , 1 ) + Data.Ek( Para.dipol_trans(ii,2),: , 1 ) ).' ;
    
    Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).' / 10; 
%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) + 1i * Data.dipol{ii}(2,:) ).';
    

%     Rausgek�rzt V
%     Bloch.feff( Bloch.ind(:,ii) ) = 1 - ( Data.fk(Para.dipol_trans(ii,1),:) + Data.fk(Para.dipol_trans(ii,2),:) ).';      
    
end

% Bloch.dipol = 5e4 * ones(Para.nr.k * Para.nr.dipol,1);                % ? 1-3 too strong.


Bloch.gamma = 10;
Bloch.E0 = 1e-7;
Bloch.t_peak = 2.038 * 1e-3;
Bloch.sigma = 1e-3;
Bloch.nrk = Para.nr.k;



% Kommt noch dazu
Emin = -2000;
Emax = 0;
E = linspace(Emin,Emax,3001)';

Bloch.w = E / constAg.hbar;             % Energiefenster in omega ???

Para.nr.w = numel(Bloch.w);



1;

%% Zeitentwicklung

Bloch.coul_ctrl = 1;                    % Coulomb Interaktion


tspan = [0 0.4];
% tspan = linspace(0, 0.4, 10000);
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
% plot( E + Data.EGap , imag(chi_w) )
% plot( E  , imag(chi_w) )
% 
% hold on
% load('spec_V_dip.mat')
% load ist_egal.mat
% plot(spec_12(:,1)*Bloch.hbar , spec_12(:,3),'r--')
% plot(spec(:,1)*Bloch.hbar , spec(:,3))


% xlim([-500 0])

Data.E = E;
Data.chi_w = chi_w;

%%
% close all
% file = 'abs_spec_0.000E+00_0.000E+00_3.000E+02_2_2_3.18_1.000E+00_1.000E+00cR_1.400E+01_60_30_1.000E-07_1.000E-03_+0.000E+00_HF_self_g0w0-tb_3_r_c_me_soc_1.519E+01_2.dat';
% Aspec = ...
%     importdata( file );
% 
figure
plot(E + Data.EGap , imag(chi_w))
% hold on
% plot(Aspec(:,1),Aspec(:,5),'r--')


%%
if Ctrl.profile_flag == 1
    profile viewer
    profile off
end
