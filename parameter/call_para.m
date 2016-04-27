function Para = call_para(Ctrl)



[Para.TB.liu.names, Para.TB.liu.values] = load_TB_parameter_old(Ctrl);

load_TB_parameter( Ctrl )





Para.symmpts{1} = {'\Gamma', 'K', 'K*', 'M'};
Para.symmpts{2} = 2 * pi / (3 * Para.TB.liu.values(1)) ...
    * [0, 0 ; 2, 0 ; 1, sqrt(3) ; 3 / 2, sqrt(3) / 2 ]';
Para.rezGV = 2 * pi / Para.TB.liu.values(1) * [1, -1 / sqrt(3); 0, 2 / sqrt(3)]';

Para.area_real = 3 * sqrt(3) / 2 * (Para.TB.liu.values(1))^2;
Para.area_BZ = 3 * sqrt(3) / 2 * (norm(Para.symmpts{2}(:,2)))^2;

Para.aBZred = norm(Para.symmpts{2}(:,2)) / Ctrl.k_mesh_mp.qr;
Para.area_sBZ = 3 * sqrt(3) / 2 * Para.aBZred^2;
Para.qmin = Para.aBZred * sqrt(3);

Para.coul_pol = 3 * sqrt(3) * log(3) * Para.aBZred / Para.area_sBZ;

Para.coul_screened = [[1.17 , 7.16 , 0.199, 2.675];
    [0.456 , 14.03 , 0.242, 2.930]; [0.46 , 13.97 , 0.24, 2.930];
    [1.288 , 6.88 , 0.232, 2.682]; [0.713 , 9.9 , 0.246, 2.855];
    [1.273 , 6.92 , 0.227, 2.682]];
Para.coul_kappa = 0;             % Kappa, because of Singularity
Para.dipol_trans = [1, 2 ; 1 , 3 ; 2, 1 ; 3 , 1 ; 4 , 5 ; 4 , 6 ; 5 , 4 ; 6 , 4 ];