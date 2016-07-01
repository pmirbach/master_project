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

% figure
% plot(k(1,:),k(2,:),'-x')
wk = pts_weight_new(k, newsymmpts, 30 * eps);
kall = red_to_BZ(k);

figure
hold on
for ii = 1:6
    plot(kall(1,:,ii),kall(2,:,ii),'x')
end


1
