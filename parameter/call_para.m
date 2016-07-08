function [ Para , W90Data ] = call_para(Ctrl, constAg)


if strcmp( Ctrl.TB_modell , 'ab_initio' ) 
    
    a0_form = num2str( 10 * Ctrl.lattice_constant , '%.3f' );
    seed = fullfile('Tight_Binding','ab_initio','02_Materials',Ctrl.material,['a0_',a0_form, 'A'],'01_WannierTB','02_G0W0','02_Mo3d','wannier90');
    SOCSettings.type  = 'none'; % 'first' or 'second' or 'none'
  
    % load wannier90 Data
    W90Data = loadW90Data(seed, SOCSettings);
    
elseif strcmp( Ctrl.TB_modell , 'liu' )
    
    [Para.TB, Para.TB_names] = load_TB_parameter( Ctrl );   % Liu parameter
    Ctrl.lattice_constant = Para.TB(1);
    W90Data = [];
    
end



Para.real.a = Ctrl.lattice_constant;
Para.real.area = 3 * sqrt(3) / 2 * ( Para.real.a )^2;



if strcmp( Ctrl.k_mesh.type , 'symm' )
    
    Para.k.GV = 2 * pi / Para.real.a * [ [ 2 / sqrt(3) ; 0 ] , [ 1 / sqrt(3) ; 1 ] ];
    
    Para.BZred.symmpts{1} = { 'K' , '\Gamma' , 'K*' , 'M' };
    Para.BZred.symmpts{2} = Para.k.GV * [ [ 1 ; 1 ] / 3 , [ 0 ; 0 ] , [ 2 ; -1 ] / 3 , [ 1 ; 0 ] / 2 ];
    
    Para.k.alpha = pi / 6;
    
    if strcmp( Ctrl.TB_modell , 'liu' )
        Para.BZred.symmpts{1} = { '\Gamma' , 'K*' , 'K' , 'M' };
    end
    
elseif strcmp( Ctrl.k_mesh.type , 'liu' )
        
    Para.k.GV = 2 * pi / Para.real.a * [ [ 1 ; -1 / sqrt(3) ] , [ 0 ; 2 / sqrt(3) ] ];
    
    Para.BZred.symmpts{1} = { '\Gamma' , 'K' , 'K*' , 'M' };
    Para.BZred.symmpts{2} = Para.k.GV * [ [ 0 ; 0 ] , [ 2 ; 1 ] / 3 , [ 1 ; 2 ] / 3 , [ 1 ; 1 ] / 2 ];
    
    Para.k.alpha = 0;
        
end



Para.energy_conversion = 1e3;                               % Umrechnung von ev in mev

Para.TB_ind{1} = [1,4]; 
Para.TB_ind{2} = [2,3,5,6];




Para.BZ.a = norm(Para.BZred.symmpts{2}(:,2));
Para.BZ.area = 3 * sqrt(3) / 2 * Para.BZ.a^2;

Para.BZsmall.a = Para.BZ.a / Ctrl.k_mesh.qr;
Para.BZsmall.area = 3 * sqrt(3) / 2 * Para.BZsmall.a^2;

Para.k.qmin = Para.BZsmall.a * sqrt(3);                                                 % Halbe???

Para.k.b = call_bnn( Para.k.GV );
Para.nr.b = size(Para.k.b,2);

[Para.coul.screened, Para.coul.names] = load_coul_parameter( Ctrl );
Para.coul.pol = 3 * sqrt(3) * log(3) * Para.BZsmall.a / Para.BZsmall.area;
Para.coul.kappa = 0;             % Kappa, because of Singularity


Para.dipol_trans = [ 1, 2 ; 4 , 5  ];
% Para.dipol_trans = [1, 2 ; 1 , 3 ; 4 , 5 ; 4 , 6 ];
Para.nr.dipol = size(Para.dipol_trans,1);

Para.coul_indices = call_coul_indices;


Para.vorf.dipol = constAg.ec;

Para.vorf.coul = constAg.ec^2 / ( 2 * constAg.eps_0 * Ctrl.eps_r);

Para.nr.tri = 6;





