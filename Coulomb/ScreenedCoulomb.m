%% Header
clear variables
clear global
close all
clc
profile on

load('KonstantenAg.mat')    % Naturkonstanten (Alex und Co)

%% Ctrl
Ctrl.material = 'MoS2';   % Materialien: MoS2 WS2 MoSe2 WSe2 MoTe2 WTe2
Ctrl.method = 'TNN';      % Möglich:   NN , TNN
Ctrl.SOC = 1;             % Spin-Orbit-Coupling

%% Material & Tight-Binding Parameter & Hochsymmetriepunkte

[Parameter.TB.liu.names, Parameter.TB.liu.values] = load_TB_parameter(Ctrl);

Parameter.symmpts{1} = {'\Gamma', 'K', 'K*', 'M'};
Parameter.symmpts{2} = 2 * pi / (3 * Parameter.TB.liu.values(1)) ...
    * [0, 0 ; 2, 0 ; 1, sqrt(3) ; 3 / 2, sqrt(3) / 2 ]';
Parameter.rezGV = 2 * pi / Parameter.TB.liu.values(1) ...
    * [1, -1 / sqrt(3); 0, 2 / sqrt(3)]';

%% Orbital Parameter

para = [[1.17 , 7.16 , 0.199, 2.675];[0.456 , 14.03 , 0.242, 2.930];
    [0.46 , 13.97 , 0.24, 2.930]; [1.288 , 6.88 , 0.232, 2.682];
    [0.713 , 9.9 , 0.246, 2.855]; [1.273 , 6.92 , 0.227, 2.682]];

%% Konstanten & (Material-) parameter
a = Parameter.TB.liu.values(1);     % Gitterkonstantek
A = 3/2 * a^2 * sqrt(3);              % Flächeninhalt der Elementarzelle

%% Vorfaktoren etc
kappa = 0.1;
vorf = constAg.ec^2 / ( 2 * constAg.eps_0 * A);
q = 0:0.001:13.2;
N_q = size(q,2);


%% Screened, all in one

figure
for jj = 1:6
    subplot(2,3,jj)
    [V,U] = fun_Coul_screened(q,para(jj,:),kappa);
    V = vorf * V;
    U = vorf * U;
    V_hartree = vorf * fun_Coul_screened_long(para(jj,:));
    V2 = V_hartree * ones(1,N_q);
    plot(q,V*1e-3)
    hold on
    plot(q,U*1e-3)
%     plot(q,V2*1e-3)
    axis([0, 13.2, 0, 30])
end

% profile report
profile off