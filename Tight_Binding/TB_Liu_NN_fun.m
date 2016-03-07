function [H_NN] = TB_Liu_NN_fun(k,Parameter)


a = Parameter.TB.liu.values(1);
eps_1 = Parameter.TB.liu.values(3);
eps_2 = Parameter.TB.liu.values(4);
t_0 = Parameter.TB.liu.values(5);
t_1 = Parameter.TB.liu.values(6);
t_2 = Parameter.TB.liu.values(7);
t_11 = Parameter.TB.liu.values(8);
t_12 = Parameter.TB.liu.values(9);
t_22 = Parameter.TB.liu.values(10);

alpha = k(1) / 2 * a;
ca = cos(alpha);
sa = sin(alpha);
beta = sqrt(3) * k(2) / 2 * a;
cb = cos(beta);
sb = sin(beta);


h_0 = 2 * t_0 * ( cos(2 * alpha) + 2 * ca * cb ) + eps_1;
h_1 = - 2 * sqrt(3) * t_2 * sa * sb ...
    + 2i * t_1 * ( sin(2 * alpha) + sa * cb );
h_2 = 2 * t_2 * ( cos(2 * alpha) - ca * cb ) ...
    + 2i * sqrt(3) * t_1 * ca * sb;
h_11 = 2 * t_11 * cos(2 * alpha) + ( t_11 + 3 * t_22 ) * ca * cb + eps_2;
h_22 = 2 * t_22 * cos(2 * alpha) + ( 3 * t_11 + t_22 ) * ca * cb + eps_2;
h_12 = sqrt(3) * ( t_22 - t_11 ) * sa * sb + 4i * t_12 * sa * ( ca - cb );

H_NN = [h_0, h_1, h_2; h_1', h_11, h_12; h_2', h_12', h_22];