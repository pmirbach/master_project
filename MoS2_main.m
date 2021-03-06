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

% System of units:
% Energy: mEv, Length: nm, Time: ps, Electric current: pA

load('KonstantenAg.mat')    % Naturkonstanten (Ag Jahnke)
[ Para , W90Data , Coul_ME ] = call_para(Ctrl, constAg);


%% k - mesh

if strcmp( Ctrl.k_mesh.type , 'symm' )
    
    [ Data.k , Data.wk , Para.nr.k , Para.k_ind ] = k_mesh_AG(Ctrl, Para);
    
elseif strcmp( Ctrl.k_mesh.type , 'liu' )
    
    [ Data.k , Data.wk , Para.nr.k , Para.k_ind.symm ] = k_mesh_liu(Ctrl, Para);
    
end


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
% [fig.bandstr_surf, fig.bandstr_path] = plot_bandstr(Ctrl,Para,Data.k,Ek_old(:,:,1)/1000,[2 3]);
% fig.bandstr_path.Position = [0.1 0.1 .4 0.6];

%% Simulation-preperations
fprintf('Preperations:              Start'); tic

[Prep.CV, Prep.CV_noSOC, Prep.minq, Prep.V_ab_interpl ] = prep( Ctrl, Para, Data, Prep.Ev_noSOC , Coul_ME );

fprintf('   -   Finished in %g seconds\n',toc)

% [fig.bandstr_surf2, fig.bandstr_path2] = plot_bandstr(Ctrl,Para,Data.k,Prep.Eks,[1 2]);


%% Dipolmatrix
fprintf('Dipol:                     Start'); tic

Data.dipol = dipol( Para , Prep , Data );

fprintf('   -   Finished in %g seconds\n',toc)


%% Dipol Plot
% titlestr = {'1 \rightarrow 2 \uparrow','1 \rightarrow 3 \uparrow','1 \rightarrow 2 \downarrow','1 \rightarrow 3 \downarrow'};
% [fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,Data.dipol,[2 2],titlestr);
% clear titlestr

A = cell2mat(Data.dipol);
As = reshape(A,Para.nr.k,[]);

% titlestr = {'1 \rightarrow 2, \uparrow, \sigma_+','1 \rightarrow 2, \downarrow, \sigma_+',...
%     '1 \rightarrow 2, \uparrow, \sigma_-','1 \rightarrow 2, \downarrow, \sigma_-'};
titlestr = {'\uparrow, \sigma_+','\downarrow, \sigma_+',...
    '\uparrow, \sigma_-','\downarrow, \sigma_-'};
[fig.dipolUp_surf, fig.dipolUp_path] = plot_dipol(Ctrl,Para,Data.k,As*1e-4,[2 2],titlestr);



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


%% Coulomb-Me for Rabi-energy renormization

fprintf('Rabi-Energie Coulombmatrixelemente:  Start'); tic

[ Data.V_rabi ] = coulomb_rabi_f_interpl(Ctrl, Para, Prep );

fprintf('   -   Finished in %g seconds\n',toc)

% %%
% % figure; 
% hold on
%%
% huhu = figure(8);
% 
% huhu.Position = [100 100 900 900];
% 
% set(gcf,'color','w')
% scatter3(Data.k(1,:,1),Data.k(2,:,1),Data.V_rabi(end,:,1)/6 , 20, (Data.V_rabi(end,:,1)/6).^(1/1) ,'filled' )
% set(gca,'fontsize',14)
% colormap jet
% 
% axis off
% hold on
% 
% cb = colorbar; 
% set(cb,'position',[.8 .2 .05 .5])
% cb.Title.String = 'eV';
% 
% corners = [Para.BZred.symmpts{2};0 0 0 0];
% plot3( corners(1,1:2),corners(2,1:2),corners(3,1:2), 'k', 'linewidth',1 )
% plot3( corners(1,2:3),corners(2,2:3),corners(3,2:3), 'k', 'linewidth',1 )
% plot3( corners(1,[3 1]),corners(2,[3 1]),corners(3,[3 1]), 'k', 'linewidth',1 )
% 
% corners_string = {'K' '\Gamma' 'K''' 'M'};
% text(corners(1,:) + 0.8 * [0.3 -0.3 0.3 0.3], ...
%         corners(2,:) + 0.5 * [-1 -1.4 1 0.5], ...
%         3000*[1 1 1 1], ...
%         corners_string,'FontSize',14)


%% Structure for all needed variables in ode

[ Bloch , Data.energy , Para.nr.w ] = call_bloch_structure( Ctrl, constAg , Para , Data );


%% Zeitentwicklung

psik_E_ini = zeros(1 , Para.nr.dipol * Para.nr.k + size(Bloch.w,1) * 2);

options=odeset('OutputFcn',@odeprog,'Events',@odeabort,'RelTol',1e-6,'AbsTol',1e-8);
ops = odeset('OutputFcn',@odetpbar,'RelTol',1e-6,'AbsTol',1e-8); 
% opts = odeset('RelTol',1e-1,'AbsTol',1e-3);
[t,psik_E] = ode113(@(t,psik_E) dgl_bloch(t,psik_E,Bloch) , Ctrl.ode.tspan , psik_E_ini , ops);
% [t,psik_E] = ode113(@(t,psik_E) dgl_bloch(t,psik_E,Bloch), tspan, psik_E_ini, options);


%%

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

Data.chi_w = chi_w;

%% Umrechnung auf Absorptionskoeffizienten

w = ( Data.energy.' + Data.EGap ) / constAg.hbar;
Y_w = 1i * w / ( 2 * constAg.c  * constAg.eps_0 ) .* Data.chi_w;
R = abs( Y_w ).^2 ./ abs( 1 - Y_w ).^2;
T = 1 ./ abs( 1 - Y_w ).^2;

alpha = 1 - R - T;

%%

figure(10)
hold on
set(gcf,'color','w')
set(gcf,'units','normalized','position',[.1 .1 .8 .45])

plot( (Data.energy )/1000 , alpha )
% plot( (Data.energy + Data.EGap)/1000 + W90Data.EGapCorr  , alpha )
set(gca,'fontsize',16)
box on

xlabel('Energie in eV')
ylabel('Absoprtionskoeff. \alpha')

%%
% xlim([min((Data.energy + Data.EGap)/1000 + W90Data.EGapCorr) max((Data.energy + Data.EGap)/1000 + W90Data.EGapCorr)])
% legend('\uparrow & \downarrow', '\uparrow', '\downarrow','location','northwest')
% linie = plot([Data.EGap/1000 + W90Data.EGapCorr , Data.EGap/1000 + W90Data.EGapCorr ],get(gca,'YLim'),'color',0.4*[1 1 1]);

%% Saving

% Ergebnis.energy = Data.energy;
% Ergebnis.alpha = alpha;
% Ergebnis.chi_w = chi_w;
% Ergebnis.EGap = Data.EGap;
% Ergebnis.EGapCorr = W90Data.EGapCorr;

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
