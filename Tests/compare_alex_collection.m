function [] = compare_alex_collection()

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

%% Vergleich Alexander #2

% Dipol_Alex = importdata('Alexander/dip.dat'); % kx,ky,12up,13up,12dwn,13dwn, alles nochmal mit anderer Pol.


%% Vergleich der Coulomb Matrixelemente mit Alex

% Coul_Alex = importdata('Alexander/Coul_Fock_60.dat'); % kx,ky,12up,13up,12dwn,13dwn, alles nochmal mit anderer Pol.
% nrkt = size(Coul_Alex,1) / 5;
% ACk = Coul_Alex(:,1:2);
% ACC = Coul_Alex(:,3);
% % 
% D = [cos(pi/6) -sin(pi/6) ; sin(pi/6) cos(pi/6)];
% 
% ACk2 = D * ACk.';
% 
% B1 = reshape( ACk2 , 2 , [] , 5 );
% ACC2 = reshape(ACC,[],5);
% 
% % color_sc = ( ACC2(:,4) / max( ACC2(:,4) ) );       
% % C = [ zeros(Para.nr.k*6,1) , color_sc , zeros(Para.nr.k*6,1) ];     % green
% % C = color_sc * [ zeros(Para.nr.k*6,1) , color_sc , zeros(Para.nr.k*6,1) ];     % green
% 
% 
% % Para.dipol_trans = [1 2];
% % Para.dipol_trans = [4 5];
% % Para.nr.dipol = size(Para.dipol_trans,1);
% 
% [Para.coul.screened, Para.coul.names] = load_coul_parameter( Ctrl );
% 
% % [V_rabi_fock] = coulomb_rabi_f_2(Ctrl, Para, Prep);  
% 
% Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);                                             % ??? Kein Unterschied???
% 
% [V_rabi_fock2] = coulomb_rabi_f_2(Ctrl, Para, Prep);  
% 
% % scatter3(Data.k(1,:),Data.k(2,:),V_rabi_fock(Para.symm_indices(2),:)-V_rabi_fock2(Para.symm_indices(2),:),10,'r+')
% 
% 
% % figure; scatter3(B1(1,:,1),B1(2,:,1),ACC2(:,2),10,'bx')
% % hold on
% % scatter3(Data.k(1,:),Data.k(2,:),V_rabi_fock2(Para.symm_indices(2),:),10,'r+')




%%

% Para.dipol_trans = [1, 2 ; 1 , 3 ; 4 , 5 ; 4 , 6 ];
% [Para.coul.screened, Para.coul.names] = load_coul_parameter( Ctrl );
% Para.coul.screened = Para.coul.screened([1,3,2,6,5,4],:);  
% [V_rabi_fock2] = coulomb_rabi_f_2(Ctrl, Para, Prep);  

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


