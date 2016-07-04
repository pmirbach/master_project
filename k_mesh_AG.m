function k_mesh_AG(Ctrl, Para)
close all

a = 0.318;

b1 = 2 * pi / a * [ 2 / sqrt(3) ; 0 ];
b2 = 2 * pi / a * [ 1 / sqrt(3) ; 1 ];


qr = Ctrl.k_mesh_mp.qr;               % Feinheit des meshes


ur = (-2*qr/3:2*qr/3) / qr;                   % Einteilung des rez. Gittervektors
A = ones(numel(ur),1) * ur;                     % Erzeugung eines meshes
k_mesh_x = A' * b1(1) + A * b2(1);
k_mesh_y = A' * b1(2) + A * b2(2);

k0 = [k_mesh_x(:), k_mesh_y(:)]';     % Alle k-Punkte


A2 = -ones(numel(ur),1) * ur;
k_mesh_x2 = A2' * b1(1) + A2 * b2(1);
k_mesh_y2 = A' * b1(2) + A * b2(2);

k02 = [k_mesh_x2(:), k_mesh_y2(:)]';     % Alle k-Punkte

kunten = k0( : , k0(2,:) < 0 );
kmitte = k0( : , k0(2,:) == 0 );
koben = k02( : , k02(2,:) > 0 );

kuntenn = kunten(: , end:-1:1);
kobenn = koben(: , end:-1:1);


k1 = [ kobenn , kmitte , kuntenn ];





% figure
% plot(k1(1,:),k1(2,:),'-x')


newsymmpts = 2 * pi / ( 3 * a ) * [ [ sqrt(3); 1 ] , [ 0; 0 ], [ sqrt(3); -1 ] ];

k = pts_triangle_fun_new(k1 , newsymmpts , 30 * eps);

Nrk = size(k,2);
Nrmitte = sum( k(2,:) == 0 );
NrSide = ( Nrk - Nrmitte ) / 2;

ind_oben = 1 : NrSide;
ind_mitte = NrSide + 1 : NrSide + Nrmitte;
ind_unten = Nrk : -1 : Nrk - NrSide + 1;


wk = pts_weight_new(k, newsymmpts, 30 * eps);



kall = zeros( [size(k),6] );

kall(:,:,1) = k;
kall(:,:,4) = kall(:,:,1);
kall(1,:,4) = - kall(1,:,4);

C3 = [ cos(2 * pi / 3), -sin(2 * pi / 3); ...
    sin(2 * pi / 3) cos(2 * pi / 3) ];

kall(:,:,3) = C3 * kall(:,:,1);
kall(:,:,5) = C3 * kall(:,:,3);
kall(:,:,6) = C3 * kall(:,:,4);
kall(:,:,2) = C3 * kall(:,:,6);





% ind = 30;
% figure; hold on
% for ii = 1:6
%     plot(kall(1,:,ii),kall(2,:,ii),'x')
%     plot(kall(1,ind,ii),kall(2,ind,ii),'ko')
% end



% % Test Plot: All in One
% 
% kind = 1:Nrk;
% 
% bool_or = false(3, Nrk);
% bool_or(1,:) = kind <= max( ind_oben );
% bool_or(2,:) = kind > max( ind_oben ) & kind < min( ind_unten );
% bool_or(3,:) = kind >= min( ind_unten );
% 
% bool_w = false(3, Nrk);
% bool_w(1,:) = wk == 6;
% bool_w(2,:) = wk == 3;
% bool_w(3,:) = wk == 1;
% 
% colors = {'b','r','g','c','m'};
% marker = {'h','x','^','+','v'};
% figure; hold on
% for ni = 1:3
%     for nj = 1:3
%         asdf = [ colors{nj}, marker{ni} ];
%         plot(k(1,bool_w(ni,:) & bool_or(nj,:),1),k(2,bool_w(ni,:) & bool_or(nj,:),1),asdf)
%     end
% end



% % Test Plot: Gewichtung
% in = wk == 6;
% on = wk == 3;
% symm = wk == 1;
% figure; hold on
% plot(kall(1,in,2),kall(2,in,2),'rx')
% plot(kall(1,on,2),kall(2,on,2),'r^')
% plot(kall(1,symm,2),kall(2,symm,2),'rs')

% % Test Plot: Indizierung 1.1 - Speicherordnung
% figure; hold on
% for ii = 1:6
%     plot(kall(1,:,ii),kall(2,:,ii),'x')
% end
% % Test Plot: Indizierung 1.2 - aequivalente Punkte
% hold on
% ind = 1;
% ind_eq = Nrk - ind + 1;
% plot(kall(1,ind,1),kall(2,ind,1),'ko')
% plot(kall(1,ind_eq,1),kall(2,ind_eq,1),'bo')

% Test Plot: Indizierung 3 - Indizes unten, mitte, oben

color_matrix = [ 160 82 45 ; 205 133 63 ; 222 184 135 ; 85 107 47 ; 107 142 35 ; 189 183 107] / 255;
opac = 1/10;

figure; hold on
for ni = 1:6
    scatter(kall(1,ind_oben,ni),kall(2,ind_oben,ni),'^','MarkerEdgeColor',color_matrix(ni,:),'MarkerFaceColor',color_matrix(ni,:),'MarkerFaceAlpha',opac)
    scatter(kall(1,ind_mitte,ni),kall(2,ind_mitte,ni),'s','MarkerEdgeColor',color_matrix(ni,:),'MarkerFaceColor',color_matrix(ni,:),'MarkerFaceAlpha',opac)
    scatter(kall(1,ind_unten,ni),kall(2,ind_unten,ni),'v','MarkerEdgeColor',color_matrix(ni,:),'MarkerFaceColor',color_matrix(ni,:),'MarkerFaceAlpha',opac)
end

testinds = [121, 250];

for ii = 1:length(testinds)
    
    for ni = 1:6
        
        indup = ind_oben( testinds( ii ) );
        inddwn = ind_unten( testinds( ii ) );
        
        plot(kall(1,indup,ni),kall(2,indup,ni),'ko')
        plot(kall(1,inddwn,ni),kall(2,inddwn,ni),'ro')
        
    end
    
end






1

% % Test Plot: Entire BZ
% figure; hold on
% for ii = 1:6
%     plot(kall(1,:,ii),kall(2,:,ii),'x')
% end


1
