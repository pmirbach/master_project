%% Header
% clear variables
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
[ Para , W90Data ] = call_para(Ctrl, constAg);

% Achtung: orbital order in Maltes TB modell different
% Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);                                             % ??? Kein Unterschied???






%% Monkhorst-Pack

[ Data.k , Data.wk , Para.nr.k , Para.symm_indices ] = k_mesh_mp(Ctrl, Para);



%%
% 
% k_mesh_AG(Ctrl, Para)


%%


file = 'Daniel/MoS2_gradtest.mat';
[Para, Data, Daniel] = load_daniel_data(Ctrl, Para, Data, file); % overwrites k, wk, coul_pol, wk_area, nr_k

% Fast kx, ky for scatter3 plots:
kx = Data.k(1,:,1);
ky = Data.k(2,:,1);




%% Tight-Binding

if strcmp(Ctrl.TB_modell,'ab_initio') 
    
    fprintf('Tight-binding (ab initio): Start'); tic 
    
    [Data.Ek, Data.Ev, Prep.Ek_noSOC, Prep.Ev_noSOC, Prep.H_grad_kx, Prep.H_grad_ky] = tight_binding_roesner(Ctrl, Para, Data, W90Data);  
    
elseif strcmp(Ctrl.TB_modell,'liu')
    
    fprintf('Tight-binding (liu):       Start'); tic 
    
    [Data.Ek, Data.Ev, Prep.Ek_noSOC, Prep.Ev_noSOC, Prep.H_grad_kx, Prep.H_grad_ky] = tight_binding_liu(Ctrl, Para, Data);
    
else
    error('TB-modell must be ab_initio or liu!')
end

fprintf('   -   Finished in %g seconds\n',toc)

[fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Data.Ek(:,:,1),[1 2]);



%% Simulation-preperations
fprintf('Preperations:              Start'); tic

[Prep.Eks, Prep.CV, Prep.CV_noSOC, Prep.minq] = prep(Para, Data, Prep.Ev_noSOC);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.bandstr_surf2, fig.bandstr_path2] = plot_bandstr(Ctrl,Para,Data.k,Prep.Eks,[1 2]);



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
% scatter3(kx,ky,Ploter.dipol_r(:,1)/10)

%% Thermische Anregung
fprintf('Excitation:                Start'); tic

[Data.fk, Para.mu] = excitation(Ctrl,constAg,Para,Data.wk,Prep.Eks);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.exc_surf, fig.exc_path] = plot_excitation(Ctrl,Para,Data.k,Data.fk,[2 3]);


%% Coulomb WW
fprintf('Coulomb matrix:            Start'); tic

[Data.V.f, Data.V.h, Data.V.h_off] = coulomb_hf( Ctrl , Para , Prep );

fprintf('   -   Finished in %g seconds\n',toc)

%% Coulomb Plots

[fig.coulomb_up, fig.coulomb_down] = plot_coulomb( Ctrl , Para,Data.k , Data.V.f , Para.symm_indices(2) );

%% Band renorm

fprintf('Band renormalization:      Start'); tic

[Data.Ek_hf, Data.Ek_h, Data.Ek_f, Test.ren_h] = renorm2(Para, Data.Ek, Data.V, Data.fk, Data.wk);

fprintf('   -   Finished in %g seconds\n',toc)

[fig.ren_bandstr_surf, fig.ren_bandstr_path] = plot_renorm_bandstr(Ctrl,Para,Data.k,[Data.Ek(:,:,1);Data.Ek_h],[2 3]);

% as1 = plot_path(Ctrl,Para,Data.k,Test.ren_h,200);      % Test der Hartree Renormierung

%%
% as2 = plot_path(Ctrl,Para,Data.k,Prep.Eks,200);

%% structure für Variablen für Blochgleichungen

% Para.dipol_trans = [1, 2 ; 1 , 3 ; 4 , 5 ; 4 , 6 ];
% Para.dipol_trans = [1 2; 4 5];
% Para.dipol_trans = [4 5];
Para.nr.dipol = size(Para.dipol_trans,1);




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
    Bloch.Esum( Bloch.ind(:,ii) ) = ( Prep.Eks( Para.dipol_trans(ii,1),: ) + Prep.Eks( Para.dipol_trans(ii,2),: ) ).' ;
%     Bloch.Esum( Bloch.ind(:,ii) ) = ( - Data.Ek( Para.dipol_trans(ii,1),: , 1 ) + Data.Ek( Para.dipol_trans(ii,2),: , 1 ) ).' ;
    
    Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) - 1i * Data.dipol{ii}(2,:) ).' / 10; 
%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) + 1i * Data.dipol{ii}(2,:) ).';
    
    Bloch.feff( Bloch.ind(:,ii) ) = 1 - ( Data.fk(Para.dipol_trans(ii,1),:) + Data.fk(Para.dipol_trans(ii,2),:) ).';
    
end

% Bloch.dipol = 5e4 * ones(Para.nr.k * Para.nr.dipol,1);                % ? 1-3 too strong.


Bloch.gamma = 10;
Bloch.E0 = 1e-7;
Bloch.t_peak = 2.038 * 1e-3;
Bloch.sigma = 1e-3;
Bloch.nrk = Para.nr.k;



% Kommt noch dazu
Emin = -1000;
Emax = 0;
E = linspace(Emin,Emax,1001)';

Bloch.w = E / constAg.hbar;             % Energiefenster in omega ???

Para.nr.w = numel(Bloch.w);

%%

% load('MoS2_631.mat')
% % load( 'MoS2_35x35.mat' ) ;
% 
% 
% 
% 
% Aint = importdata('test_1.dat');
% 
% Akx = Aint(:,1);
% Aky = Aint(:,2);
% % 
% 
% k1 = [ Akx , Aky ];
% D = [cos(pi/6) -sin(pi/6) ; sin(pi/6) cos(pi/6)];
% 
% k2 = D * k1.';
% 
% 
% 
% Pint = 1 / ( 2 * pi )^2 * Bloch.coulomb(:,:,1) * ( Bloch.wkentire );
% Pint2 = 1 / ( 2 * pi )^2 * Bloch.coulomb(:,:,2) * ( Bloch.wkentire );
% 
% % scatter3( MoS2.kpts(:,1,1) , MoS2.kpts(:,2,1) , sum( squeeze( MoS2.COUL.D(1,:,:) ) , 2 )  )
% % hold on
% % scatter3(kx,ky,Pint )
% 
% 
% % scatter3(kx,ky,Pint-  sum( squeeze( MoS2.COUL.D(1,:,:) ) , 2 ))
% % 
% % 
% % plot(kx(1:10),ky(1:10),'bx-')
% % hold on
% % plot(MoS2.kpts(1:10,1,1),MoS2.kpts(1:10,2,1),'rx-')
% 
% A = sum( squeeze( MoS2.COUL.D(1,:,:) ) , 2 );
% B = sum( squeeze( MoS2.COUL.D(3,:,:) ) , 2 );
% 
% 
% Pkall = round( [ Data.k(1,:,1) ; Data.k(2,:,1) ] , 7 );
% Dkall = round(MoS2.kpts(:,1:2,1),7);
% % figure
% % compare_alex( Pkall, Dkall , Pint, A, 'abs' )
% % figure
% % compare_alex( Pkall, Dkall , Pint2, B, 'abs' )
% 
% 
% 
% compare_alex( Pkall, Dkall , Ploter.dipol_r(:,1)/10, abs(MoS2.d_k(3,:)), 'abs' )
% % 
% % figure; scatter3( Dkall(:,1,1), Dkall(:,2,1), abs(MoS2.d_k(1,:)) )
% % figure; scatter3(kx,ky,Ploter.dipol_r(:,1)/10)
% 
% % % close all
% % % figure
% % % scatter3( k2(1,:), k2(2,:), Aint(:,3)  )
% % % hold on
% % % scatter3(kx,ky,Pint)
% % 
% % Pkall = round( [ Data.k(1,:,1) ; Data.k(2,:,1) ] , 7 );
% % Akall = round( k2 , 7 );
% % 
% % compare_alex( Pkall, Akall , Pint, Aint(:,3), 'rel' )
% 
% 
% % Dipol_Alex = importdata('Alexander/dip_60.dat'); % kx,ky,12up,13up,12dwn,13dwn, alles nochmal mit anderer Pol.
% % Akx = Dipol_Alex(:,1);
% % Aky = Dipol_Alex(:,2);
% % 
% % k0 = [Akx.'; Aky.'];
% % 
% % k1 = [ Akx , Aky ];
% % D = [cos(pi/6) -sin(pi/6) ; sin(pi/6) cos(pi/6)];
% % 
% % k2 = D * k1.';
% % 
% % close all
% % scatter3(kx,ky,Ploter.dipol_r(:,1)/10, 'b')
% % hold on
% % scatter3(k2(1,:),k2(2,:),Dipol_Alex(:,7) , 'g')
% % scatter3( MoS2.kpts(:,1,1) , MoS2.kpts(:,2,1) , abs(MoS2.d_k(3,:)) , 'r' )




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
plot(E , imag(chi_w))

hold on
% load('spec_V_dip.mat')
% load ist_egal.mat
% plot(spec_12(:,1)*Bloch.hbar , spec_12(:,3),'r--')
% plot(spec(:,1)*Bloch.hbar , spec(:,3))


% xlim([-500 0])

%%
% % close all
% file = 'abs_spec_0.000E+00_0.000E+00_3.000E+02_2_2_3.18_1.000E+00_1.000E+00cR_1.400E+01_60_30_1.000E-07_1.000E-03_+0.000E+00_HF_self_g0w0-tb_3_r_c_me_soc_1.519E+01_2.dat';
% Aspec = ...
%     importdata( file );
% 
% figure
% plot(E + 2640.47 , imag(chi_w))
% hold on
% plot(Aspec(:,1),Aspec(:,5))


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
