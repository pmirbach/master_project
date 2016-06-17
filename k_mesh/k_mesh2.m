function k_mesh2( Ctrl, Para )

qr = Ctrl.k_mesh_mp.qr;             % Feinheit des meshes

b1 = abs( Para.k.GV(:,1) );          % Reziproke Gittervektoren
b2 = Para.k.GV(:,2);                 % Maltes G



n1 = qr * 2 / 3 ;
n2 = qr / 3 ;


r_b1 = 0 : n1 ;
r_b2 = -n2 : n2 ;

[M_b1, M_b2] = meshgrid(r_b1,r_b2);



asdf = [M_b1(:) , M_b2(:) ]';


borders = [ 0 , n1 , n2 ; 0 , -n2 , n2 ];

k = pts_triangle_fun( asdf , borders , 30 * eps );

wk = pts_weight( k , borders , 30 * eps );

1
11



% % Complete parallelogramm:
% a1 = 0:20;
% a2 = -10:10;
% [A1, A2] = meshgrid(a1,a2);
% kx = A1 / 30 * b1(1) + A2 / 30 * b2(1);
% ky = A1 / 30 * b1(2) + A2 / 30 * b2(2);
% plot(kx(:),ky(:),'x')