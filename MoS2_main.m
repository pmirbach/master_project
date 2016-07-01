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
[ Para , W90Data ] = call_para(Ctrl, constAg);

% Achtung: orbital order in Maltes TB modell different
% Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);                                             % ??? Kein Unterschied???






%% Monkhorst-Pack

[ Data.k , Data.wk , Para.nr.k , Para.symm_indices ] = k_mesh_mp(Ctrl, Para);


% k_mesh_AG(Ctrl, Para)


file = 'Daniel/MoS2_gradtest.mat';
[Para, Data, Daniel] = load_daniel_data(Ctrl, Para, Data, file); % overwrites k, wk, coul_pol, wk_area, nr_k

% Fast kx, ky for scatter3 plots:
kx = Data.k(1,:,1);
ky = Data.k(2,:,1);

%% Vergleich Alexander #1
% close all

Dipol_Alex = importdata('Alexander/dip_60.dat'); % kx,ky,12up,13up,12dwn,13dwn, alles nochmal mit anderer Pol.
Akx = Dipol_Alex(:,1);
Aky = Dipol_Alex(:,2);
% 
k0 = [Akx.'; Aky.'];

k1 = [ Akx , Aky ];
D = [cos(pi/6) -sin(pi/6) ; sin(pi/6) cos(pi/6)];

k2 = D * k1.';

% plot(kx,ky,'bx')
% hold on
% plot(k2(1,:),k2(2,:),'rx')


% Lokalisierung der Hochsymmetriepunkte bei Alex:
% K+ = (11.4075; -6.5861), K- = (11.4075; 6.5861)
kf = round(k0,4);

G = repmat([0; 0],1,size(kf,2));
Kp = repmat([11.4075; -6.5861],1,size(kf,2));
Km = repmat([11.4075;  6.5861],1,size(kf,2));

Gind = find(all(kf == G,1));
Kmind = find(all(kf == Km,1));
Kpind = find(all(kf == Kp,1));

% plot( k2(1,Gind),k2(2,Gind), 'ko' )
% plot( k2(1,Kmind),k2(2,Kmind), 'k^' )
% plot( k2(1,Kpind),k2(2,Kpind), 'ks' )


% figure;
% set(gcf,'color','w');
% plot(Akx,Aky,'-x','Color', 0 * [1 1 1])
% axis equal
% % hold on
% % plot(k2(1,:),k2(2,:),'Color', 0.8 * [0 1 1])
% grid off
% axis off




% Data.k = zeros(size(Data.k));
% Data.k(:,:,1) = k0;

% symmneu = Para.BZred.symmpts{2}(:,[1,3]);
% symmneu = [symmneu , [ -symmneu(1,2); symmneu(2,2) ]];
% wk = pts_weight(k0, symmneu, 30 * eps);
% 
% k = red_to_BZ(k0);

Energy_Alex = importdata('Alexander/disp_120.dat');

Ak0 = Energy_Alex(:,1:2);

Ak2 = D * Ak0.';




% % Vergleich Dispersion mit qr = 120 !
% figure
% scatter3(kx,ky,Prep.Eks(1,:))
% hold on
% scatter3(Ak2(1,:),Ak2(2,:),Energy_Alex(:,3))




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


%%
% close all

% % Prep.Eks(Para.TB_ind{2},:) = Prep.Eks(Para.TB_ind{2},:) - 0.1516 * 1e3;
% % Eks_corr = Prep.Eks;
% % Eks_corr(Para.TB_ind{2},:) = Eks_corr(Para.TB_ind{2},:) - 0.1516 * 1e3;
% 
% index = [3 5 7 4 6 8];
% figure
% for ii = 1:6
%     subplot(2,3,ii)
%     
%     compare_alex( Data.k(:,:,1), Ak2 , Prep.Eks(ii,:), Energy_Alex(:,index(ii)) , 'rel' )
% end

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


%% Vergleich mit Alexander
% % Vergleich Dipol Matrixelemente mit qr = 60 oder qr = 120
% figure
% for ii = 1:4
%     subplot(2,2,ii)
%     scatter3(kx,ky,Ploter.dipol_l(:,ii)/10)
%     hold on
%     scatter3(k2(1,:),k2(2,:),Dipol_Alex(:,ii+2))
% end
% 
% 
% figure
% for ii = 1:4
%     subplot(2,2,ii)
%     scatter3(kx,ky,Ploter.dipol_r(:,ii)/10)
%     hold on
%     scatter3(k2(1,:),k2(2,:),Dipol_Alex(:,ii+6))
% end

% figure; compare_alex( Data.k(:,:,1), k2 , Ploter.dipol_l(:,1)/10, Dipol_Alex(:,1+2) )

% figure
% for ii = 1:4
%     subplot(2,2,ii)
%     compare_alex( Data.k(:,:,1), k2 , Ploter.dipol_l(:,ii)/10, Dipol_Alex(:,ii+2),'rel' )
% end
% 
% figure
% for ii = 1:4
%     subplot(2,2,ii)
%     compare_alex( Data.k(:,:,1), k2 , Ploter.dipol_r(:,ii)/10, Dipol_Alex(:,ii+6),'rel' )
% end


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



%% Vergleich Alexander #2

% Dipol_Alex = importdata('Alexander/dip.dat'); % kx,ky,12up,13up,12dwn,13dwn, alles nochmal mit anderer Pol.


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
    
%     Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{dipolnr}(1,:) - 1i * Data.dipol{dipolnr}(2,:) ).';
    Bloch.dipol( Bloch.ind(:,ii) ) = 1 / sqrt(2) * abs( Data.dipol{ii}(1,:) + 1i * Data.dipol{ii}(2,:) ).' / 10; 
%     Bloch.dipol( (ii-1) * Para.nr.k + 1 : ii * Para.nr.k ) = abs( Data.dipol{ii}(1,:) ).'; 
    
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
Emax = 200;
E = linspace(Emin,Emax,1001)';

Bloch.w = E / constAg.hbar;             % Energiefenster in omega ???

Para.nr.w = numel(Bloch.w);





%% Vergleich der Coulomb Matrixelemente mit Alex

Coul_Alex = importdata('Alexander/Coul_Fock_60.dat'); % kx,ky,12up,13up,12dwn,13dwn, alles nochmal mit anderer Pol.
nrkt = size(Coul_Alex,1) / 5;
ACk = Coul_Alex(:,1:2);
ACC = Coul_Alex(:,3);
% 
D = [cos(pi/6) -sin(pi/6) ; sin(pi/6) cos(pi/6)];

ACk2 = D * ACk.';

B1 = reshape( ACk2 , 2 , [] , 5 );
ACC2 = reshape(ACC,[],5);

% color_sc = ( ACC2(:,4) / max( ACC2(:,4) ) );       
% C = [ zeros(Para.nr.k*6,1) , color_sc , zeros(Para.nr.k*6,1) ];     % green
% C = color_sc * [ zeros(Para.nr.k*6,1) , color_sc , zeros(Para.nr.k*6,1) ];     % green


% Para.dipol_trans = [1 2];
% Para.dipol_trans = [4 5];
% Para.nr.dipol = size(Para.dipol_trans,1);

[Para.coul.screened, Para.coul.names] = load_coul_parameter( Ctrl );

% [V_rabi_fock] = coulomb_rabi_f_2(Ctrl, Para, Prep);  

Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);                                             % ??? Kein Unterschied???

[V_rabi_fock2] = coulomb_rabi_f_2(Ctrl, Para, Prep);  

% scatter3(Data.k(1,:),Data.k(2,:),V_rabi_fock(Para.symm_indices(2),:)-V_rabi_fock2(Para.symm_indices(2),:),10,'r+')


% figure; scatter3(B1(1,:,1),B1(2,:,1),ACC2(:,2),10,'bx')
% hold on
% scatter3(Data.k(1,:),Data.k(2,:),V_rabi_fock2(Para.symm_indices(2),:),10,'r+')




%%

Para.dipol_trans = [1, 2 ; 1 , 3 ; 4 , 5 ; 4 , 6 ];
[Para.coul.screened, Para.coul.names] = load_coul_parameter( Ctrl );
Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);  
[V_rabi_fock2] = coulomb_rabi_f_2(Ctrl, Para, Prep);  

%%
% close all
% 
% Pkall = [Data.k(1,:);Data.k(2,:)];
% Akall = [B1(1,:,1);B1(2,:,1)];
% 
% figure
% subplot(1,2,1)
% PVw = V_rabi_fock2(Para.symm_indices(2),:,:,1);
% PV = PVw(:);
% AV = ACC2(:,2);
% compare_alex( Pkall, Akall , PV, AV, 'abs' )
% 
% subplot(1,2,2)
% PVw = V_rabi_fock2(Para.symm_indices(2),:,:,2);
% PV = PVw(:);
% AV = ACC2(:,4);
% compare_alex( Pkall, Akall , PV, AV, 'abs' )


%%


% [Pkalls, Pindex] = sortrows(Pkall);
% [Akalls, Aindex] = sortrows(Akall);
% 
% % any( any( Pkalls - Akalls ) )
% 
% 
% PV = V_rabi_fock2(Para.symm_indices(2),:);
% AV = ACC2(:,2);
% 
% PVs = PV(Pindex);
% AVs = AV(Aindex);
% 
% 
% % figure; scatter3( Akalls(:,1),Akalls(:,2),AVs,10,'bx' )
% % hold on
% % scatter3( Pkalls(:,1),Pkalls(:,2),PVs,10,'r+' )
% 
% 
% figure
% scatter3( Pkalls(:,1),Pkalls(:,2),PVs-AVs.',10,'r+' )
% 
% % plot( Data.k(1,:),Data.k(2,:),'rx' )
% 
% 
% 
% 
% 
% 
% 
% % figure
% % hold on
% % for ii = 1:6
% %     plot(B1(1,ii:6:end,1),B1(2,ii:6:end,1),'x')
% % end
% 
% 
% % ACkn = reshape( ACk' , 2 , [] , 6 , 5 );
% % plot(ACkn(1,:,1,1),ACkn(2,:,1,1),'rx')
% 
% % ACkt(:,:,1) = [ACk(1:6:end,1),ACk(1:6:end,2)];
% 
% 
% % plot(ACknt(1,:,1),ACknt(2,:,1),'rx')
% 
% % plot(ACk(1:6:end,1),ACk(1:6:end,2),'rx')
% % plot(ACk2(1,:),ACk2(2,:),'rx')
% 
% % figure; scatter3(ACk(:,1),ACk(:,2),Coul_Alex(:,3))
% % figure; scatter3(ACk2(1,:),ACk2(2,:),Coul_Alex(:,3))
% 
% 
% % figure; scatter3(ACk(1:nrkt,1),ACk(1:nrkt,2),Coul_Alex(1:nrkt,3))

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
close all
file = 'abs_spec_0.000E+00_0.000E+00_3.000E+02_2_2_3.18_1.000E+00_1.000E+00cR_1.400E+01_60_30_1.000E-07_1.000E-03_+0.000E+00_HF_self_g0w0-tb_3_r_c_me_soc_1.519E+01_2.dat';
Aspec = ...
    importdata(file);

figure
plot(E + 2640.47 , imag(chi_w))
hold on
plot(Aspec(:,1),Aspec(:,5))


%%
Aint = importdata('test_1.dat');

Akx = Aint(:,1);
Aky = Aint(:,2);
% 

k1 = [ Akx , Aky ];
D = [cos(pi/6) -sin(pi/6) ; sin(pi/6) cos(pi/6)];

k2 = D * k1.';



Pint = 1 / ( 2 * pi )^2 * Bloch.coulomb(:,:,1) * ( Bloch.wkentire );

close all
figure
scatter3( k2(1,:), k2(2,:), Aint(:,3)  )
hold on
scatter3(kx,ky,Pint)


Pint_n = zeros(Para.nr.k,1);
for nk = 1:Para.nr.k
    
    for nks = 1:Para.nr.k
        
        for tri = 1:6
        
            Pint_n(nk) = Pint_n(nk) + 1 / ( 2 * pi )^2 * Data.wk(nks) / 6 * Para.BZsmall.area * V_rabi_fock2(nk,nks,tri,1) * Para.vorf.coul;
        
        end
        
    end
    
end

figure
scatter3( k2(1,:), k2(2,:), Aint(:,3)  )
hold on
scatter3(kx,ky,Pint_n)





%%

RenKs = 0;
Ev = Prep.Ev_noSOC;
para_map = [1 2 3 ; 2 4 5 ; 3 5 6];

k0x = Data.k(1,end,1);
k0y = Data.k(2,end,1);

for nks = 1:Para.nr.k
    for tri = 1:6
        
        ksx = Data.k(1,nks,tri);
        ksy = Data.k(2,nks,tri);
        
        qb = zeros(1,Para.nr.b);
        
        for nb = 1:Para.nr.b
            
            qb(nb) = sqrt( ( k0x - ksx + Para.k.b(1,nb) )^2 + ( k0y - ksy + Para.k.b(1,nb) )^2 );
            
        end
        
        q = min( qb );
        
        if q < 1e-4
            q = 0;
        end
        
        
        Coul = 0;
        for a = 1:3
            for b = 1:3
                
                V = final_coul_scr(Prep.minq(end,nks,tri),Para.coul.screened(para_map(a,b),:),Para.coul.pol);
                               
                Coul = conj( Ev(a,2,end,1) ) * conj( Ev(b,1,nks,tri) ) * Ev(b,1,end,1) * Ev(a,2,nks,tri) * V;
                
            end
        end
        
        Coul = Para.vorf.coul * abs( Coul );
        
        RenKs = RenKs + 1 / ( 2 * pi )^2 * Data.wk(nks) / 6 * Para.BZsmall.area * Coul;
        
    end
end




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
if Ctrl.profile_flag == 1
    profile viewer
    profile off
end
