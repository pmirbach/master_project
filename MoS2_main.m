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

if strcmp( Ctrl.k_mesh.type , 'symm' )
    
    [ Data.k , Data.wk , Para.nr.k , Para.k_ind ] = k_mesh_AG(Ctrl, Para);
    
elseif strcmp( Ctrl.k_mesh.type , 'liu' )
    
    [ Data.k , Data.wk , Para.nr.k , Para.k_ind.symm ] = k_mesh_liu(Ctrl, Para);
    
end


%%




%%


% file = 'Daniel/MoS2_gradtest.mat';
% [Para, Data, Daniel] = load_daniel_data(Ctrl, Para, Data, file); % overwrites k, wk, coul_pol, wk_area, nr_k
% 
% % Fast kx, ky for scatter3 plots:
kx = Data.k(1,:,1);
ky = Data.k(2,:,1);




%% Tight-Binding

if strcmp(Ctrl.TB_modell,'ab_initio') 
    
    fprintf('Tight-binding (ab initio): Start'); tic 
    
    [Data.Ek, Data.Ev, Data.EGap, Prep.Ev_noSOC, Prep.H_grad_kx, Prep.H_grad_ky, Ek_old] = tight_binding_roesner(Ctrl, Para, Data.k, W90Data);  
    
elseif strcmp(Ctrl.TB_modell,'liu')
    
    fprintf('Tight-binding (liu):       Start'); tic 
    
    [Data.Ek, Data.Ev, Data.EGap, Prep.Ek_noSOC, Prep.Ev_noSOC, Prep.H_grad_kx, Prep.H_grad_ky] = tight_binding_liu(Ctrl, Para, Data.k);
    
else
    error('TB-modell must be ab_initio or liu!')
end

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

% titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow'};
% [fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol,[2 2],titlestr);
% clear titlestr

Ploter.dipol = zeros( Para.nr.k, Para.nr.dipol );
for ii = 1:Para.nr.dipol
    Ploter.dipol_l(:,ii) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) + 1i * Data.dipol{ii}(2,:) ).'; 
    Ploter.dipol_r(:,ii) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).'; 
end

%%

test = [Ploter.dipol_l,Ploter.dipol_r];
titlestr = {'1 \rightarrow 2, \uparrow, \sigma_+','1 \rightarrow 2, \downarrow, \sigma_+',...
    '1 \rightarrow 2, \uparrow, \sigma_-','1 \rightarrow 2, \downarrow, \sigma_-'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,test,[2 2],titlestr);

1

%%


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
Bloch.ind = reshape( 1 : Para.nr.dipol * Para.nr.k ,[], Para.nr.dipol);

% tic
% [V_rabi_fock] = coulomb_rabi_f(Ctrl, Para, Prep, Data.Ev);                           % All coulomb matrices. (in 3rd dimension)
% toc

fprintf('Rabi-Energie Coulombmatrixelemente:  Start'); tic

[V_rabi_fock_interpl] = coulomb_rabi_f_interpl(Ctrl, Para, Prep );

fprintf('   -   Finished in %g seconds\n',toc)

% %%
% % figure; 
% hold on
% scatter3(kx,ky,V_rabi_fock_interpl(1,:,1) )
% scatter3(kx,ky,V_rabi_fock(1,:,1),'r')
% % % figure; scatter3(kx,ky,V_rabi_fock_interpl(1,:,1) )

%%
% figure(15)
% subplot(1,2,1)
% plot(kx,ky,'x')
% 
% load k_stdmp_631.mat
% D = [cos(pi/6) sin(pi/6) ; -sin(pi/6) cos(pi/6)];
% k_daniel_rot = k(:,1:2,1);
% k_daniel = D * k_daniel_rot.';
% 
% subplot(1,2,2)
% plot(k_daniel(1,:),k_daniel(2,:),'x')

%%
% 
% load COUL_malwieder.mat
% 
% asdf = squeeze(COUL.D(1,end,:)) ./ squeeze(k(:,3,1));
% 
% figure(16)
% scatter3(k_daniel(1,:),k_daniel(2,:),asdf)
% hold on
% scatter3(kx,ky,V_rabi_fock_interpl(1,:,1) / 6 )
% 
% %%
% 
% figure(17)
% compare_alex( Data.k(:,:,1), k_daniel, V_rabi_fock_interpl(1,:,1) / 6 , asdf , 'abs' )

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
    
%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).' * 10; 
    
    Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).' / 10; 



%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) + 1i * Data.dipol{ii}(2,:) ).';
    

%     Rausgek�rzt V
%     Bloch.feff( Bloch.ind(:,ii) ) = 1 - ( Data.fk(Para.dipol_trans(ii,1),:) + Data.fk(Para.dipol_trans(ii,2),:) ).';      
    
end

% Bloch.dipol = 5e4 * ones(Para.nr.k * Para.nr.dipol,1);                % ? 1-3 too strong.


Bloch.gamma = 10;           % Dephrasing


Bloch.E0 = 1e-7;
Bloch.t_peak = 2.038 * 1e-3;
Bloch.sigma = 1e-3;
Bloch.nrk = Para.nr.k;



% Kommt noch dazu
Emin = -1000;
Emax = 0;
E = linspace(Emin,Emax,3001)';

Bloch.w = E / constAg.hbar;             % Energiefenster in omega ???

Para.nr.w = numel(Bloch.w);



1;

%%

kx;
ky;



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

figure
% plot( E + Data.EGap , imag(chi_w) )
plot( E  , imag(chi_w) )
% 
% hold on
% load('spec_V_dip.mat')
% load ist_egal.mat
% plot(spec_12(:,1)*Bloch.hbar , spec_12(:,3),'r--')
% plot(spec(:,1)*Bloch.hbar , spec(:,3))


% xlim([-500 0])

Data.E = E;
Data.chi_w = chi_w;

%% Umrechnung auf Absorptionskoeffizienten

w = ( Data.E.' + Data.EGap ) / constAg.hbar;
Y_w = 1i * w / ( 2 * constAg.c  * constAg.eps_0 ) .* Data.chi_w;
R = abs( Y_w ).^2 ./ abs( 1 - Y_w ).^2;
T = 1 ./ abs( 1 - Y_w ).^2;

alpha = 1 - R - T;

plot( E  , alpha )

% hold on
% filename_spec = ['abs_spec_0.000E+00_0.000E+00_3.000E+02_2_2_3.18_' ...
%     '1.000E+00_1.000E+00cR_60_30_1.000E-07_1.000E-03_+0.000E+00_HF_g0w0-tb_3_r_c_me_sock_7.596E+00_5.dat'];
% absspec0 = import_alex_spec(filename_spec);
% plot(absspec0(:,1),absspec0(:,5))

%% Saving

Ergebnis.E = E;
Ergebnis.alpha = alpha;
Ergebnis.chi_w = chi_w;
Ergebnis.EGap = Data.EGap;


%%
% plot(E + Data.EGap , imag(chi_w))

% 
% 
% % figure
% % set(gcf,'color','w');
% hold on
% 
% plot(E , imag(chi_w))
% set(gca,'fontsize',18)
% 
% xlabel('Energie E-E_{Gap} in meV')
% ylabel('\Im{\chi(\omega)}')
% % legend('Spin \uparrow','Spin \downarrow')

% title('WS_2')


%%
% warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
% measuredperaks = findpeaksG(E,imag(chi_w),.0001,2,27,18,3);
% warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
% 
% fprintf('Anzahl k-Punkte: %.0f\n',Para.nr.k)
% fprintf('A-Exziton:       %.1f\n',measuredperaks(1,2))
% fprintf('B-Exziton:       %.1f\n',measuredperaks(2,2))
% 
% 
% fileID = fopen('Ergebnisse\Konvergenz_WSe2.txt','a');
% fprintf(fileID,'%.0f %.0f %.1f %.1f\n',Ctrl.k_mesh.qr,Para.nr.k,measuredperaks(1,2),measuredperaks(2,2));
% fclose(fileID);


%%
% hold on
% file = 'abs_spec_0.000E+00_0.000E+00_3.000E+02_2_2_3.18_1.000E+00_1.000E+00cR_1.400E+01_60_30_1.000E-07_1.000E-03_+0.000E+00_HF_self_g0w0-tb_3_r_c_me_soc_1.519E+01_2.dat';
% Aspec = ...
%     importdata( file );

% hold on
% plot(Aspec(:,1),Aspec(:,5),'r--')

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
% 
% % E_t2 = Bloch.E0 * exp(-1/2 * ( ( t.' - Bloch.t_peak ) / Bloch.sigma ).^2 * 4 * log(2) );
% %
% % P_w2 = fft(P_t2);
% % E_w2 = fft(E_t2);
% %
% % chi_w2 = P_w2 ./ E_w2;
% %
% % close all
% % plot(imag(P_w2))

%%

% figure
% plot(t,real(P_t2))
% hold on
% plot(t,imag(P_t2),'r')
% legend('real','imag')

%%
% plot_surf(Para,Data.k(:,:,1),Bloch.dipol,[1,1])

%%
if Ctrl.profile_flag == 1
    profile viewer
    profile off
end
