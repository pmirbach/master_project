function Para = call_para(Ctrl, constAg)


[Para.TB, Para.TB_names] = load_TB_parameter( Ctrl );

Para.TB_ind{1} = [1,4]; 
Para.TB_ind{2} = [2,3,5,6];

Para.real.a = Para.TB(1);
Para.real.area = 3 * sqrt(3) / 2 * ( Para.real.a )^2;

Para.BZred.symmpts{1} = {'\Gamma', 'K', 'K*', 'M'};
Para.BZred.symmpts{2} = 2 * pi / ( 3 * Para.real.a ) * [0, 0 ; 2, 0 ; 1, sqrt(3) ; 3 / 2, sqrt(3) / 2 ]';

Para.BZ.a = norm(Para.BZred.symmpts{2}(:,2));
Para.BZ.area = 3 * sqrt(3) / 2 * Para.BZ.a^2;

Para.BZsmall.a = Para.BZ.a / Ctrl.k_mesh_mp.qr;
Para.BZsmall.area = 3 * sqrt(3) / 2 * Para.BZsmall.a^2;

Para.k.GV = 2 * pi / Para.real.a * [1, -1 / sqrt(3); 0, 2 / sqrt(3)]';
Para.k.qmin = Para.BZsmall.a * sqrt(3);

Para.k.b = call_bnn( Para.k.GV );
Para.nr.b = size(Para.k.b,2);

[Para.coul.screened, Para.coul.names] = load_coul_parameter( Ctrl );
Para.coul.pol = 3 * sqrt(3) * log(3) * Para.BZsmall.a / Para.BZsmall.area;
Para.coul.kappa = 0;             % Kappa, because of Singularity


Para.dipol_trans = [1, 2 ; 1 , 3 ; 4 , 5 ; 4 , 6 ];
Para.coul_indices = call_coul_indices;


Para.vorf.dipol = constAg.ec;
Para.vorf.coul = constAg.ec^2 / ( 2 * constAg.eps_0);





