function [] = get_background_screening( Para , Coul_ME , eps_2 , eps_3 , maxq )


Resta_fit_eps = 1;



minq = 0.1;           % Pol? Knapp hinter 0 anfangen?
del = 1000;

q = minq : (maxq - minq)/del : maxq;
nr_q = numel( q );


A = Para.real.area;



U_diag = zeros(3,3,nr_q);

U_diag(1,1,:) = 3 * 1 ./ q .* 1 ./ ( 1 + Coul_ME.U.gamma .* q + Coul_ME.U.delta .* q.^2 );
U_diag(2,2,:) = repmat( Coul_ME.U.micro(1) / Para.vorf.coul * A , 1 , nr_q );
U_diag(3,3,:) = repmat( Coul_ME.U.micro(2) / Para.vorf.coul * A , 1 , nr_q );


if Resta_fit_eps
    
    a = Coul_ME.eps.Resta(1);
    b = Coul_ME.eps.Resta(2);
    c = Coul_ME.eps.Resta(3);
    d = Coul_ME.eps.Resta(4);
    e = Coul_ME.eps.Resta(5);
    
    eps_inf_q = ( a + q.^2 ) ./ ( a / b  * sin( q * c ) ./ ( q * c ) + q.^2 ) + e;
    height = d;
    
else
    
    eps_inf_q = repmat( Coul_ME.eps.inf , 1 , nr_q );
    height = Coul_ME.eps.L_2d_macro;
    
end

g_2 = ( eps_inf_q - eps_2 ) ./ ( eps_inf_q + eps_2 );
g_3 = ( eps_inf_q - eps_3 ) ./ ( eps_inf_q + eps_3 );


eps_diag = zeros(3,3,nr_q);

eps_diag(1,1,:) = eps_inf_q .* ( 1 - g_2 .* g_3 .* exp( -2 * q * height ) ) ...
    ./ ( 1 + ( g_2 + g_3 ) .* exp( -q * height ) + g_2 .* g_3 .* exp( -2 * q * height ) );
eps_diag(2,2,:) = repmat( Coul_ME.eps.micro(1) , 1 , nr_q );
eps_diag(3,3,:) = repmat( Coul_ME.eps.micro(2) , 1 , nr_q );



V_ab_q = zeros(3,3,nr_q);

for ii = 1:nr_q
    
    V_diag = U_diag(:,:,ii) / eps_diag(:,:,ii);
    
    V_ab_q(:,:,ii) = Coul_ME.U_Ev * V_diag / Coul_ME.U_Ev ;
    
end

V_new = reshape( V_ab_q, 9 , [] );

figure
for ii = 1:9
    subplot(3,3,ii)
    plot(q, V_new(ii,:) )
end


1