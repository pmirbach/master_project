function [H_TNN] = TB_Liu_TNN_fun(k,Parameter)


a = Parameter(1);
eps_1 = Parameter(3);
eps_2 = Parameter(4);
t_0 = Parameter(5);
t_1 = Parameter(6);
t_2 = Parameter(7);
t_11 = Parameter(8);
t_12 = Parameter(9);
t_22 = Parameter(10);
r_0 = Parameter(11);
r_1 = Parameter(12);
r_2 = Parameter(13);
r_11 = Parameter(14);
r_12 = Parameter(15);
u_0 = Parameter(16);
u_1 = Parameter(17);
u_2 = Parameter(18);
u_11 = Parameter(19);
u_12 = Parameter(20);
u_22 = Parameter(21);



alpha = k(1) / 2 * a;
ca = cos(alpha);
c2a = cos(2 * alpha);
c3a = cos(3 * alpha);
c4a = cos(4 * alpha);
sa = sin(alpha);
s2a = sin(2 * alpha);
s3a = sin(3 * alpha);

b = sqrt(3) * k(2) / 2 * a;
cb = cos(b);
c2b = cos(2 * b);
sb = sin(b);
s2b = sin(2 * b);

v_0 = eps_1 + 2 * t_0 * ( 2 * ca * cb + c2a ) ...
    + 2 * r_0 * ( 2 * c3a * cb + c2b ) ...
    + 2 * u_0 * ( 2 * c2a * c2b + c4a );

real_v_1 = - 2 * sqrt(3) * t_2 * sa * sb ...
    + 2 * ( r_1 + r_2 ) * s3a * sb ...
    - 2 * sqrt(3) * u_2 * s2a * s2b;
imag_v_1 = 2 * t_1 * sa * ( 2 * ca + cb ) ...
    + 2 * (r_1 - r_2) * s3a * cb ...
    + 2 * u_1 * s2a * ( 2 * c2a + c2b );
v_1 = real_v_1 + 1i * imag_v_1;

real_v_2 = 2 * t_2 * ( c2a - ca * cb ) ...
    - 2 / sqrt(3) * ( r_1 + r_2 ) * ( c3a * cb - c2b ) ...
    + 2 * u_2 * ( c4a - c2a * c2b );
imag_v_2 = 2 * sqrt(3) * t_1 * ca * sb ...
    + 2 / sqrt(3) * sb * ( r_1 - r_2 ) * ( c3a + 2 * cb ) ...
    + 2 * sqrt(3) * u_1 * c2a * s2b;
v_2 = real_v_2 + 1i * imag_v_2;

v_11 = eps_2 + ( t_11 + 3 * t_22 ) * ca * cb + 2 * t_11 * c2a ...
    + 4 * r_11 * c3a * cb + 2 * ( r_11 + sqrt(3) * r_12 ) * c2b ...
    + ( u_11 + 3 * u_22 ) * c2a * c2b + 2 * u_11 * c4a; 

real_v_12 = sqrt(3) * ( t_22 - t_11 ) * sa * sb + 4 * r_12 * s3a * sb ...
    + sqrt(3) * (u_22 - u_11) * s2a * s2b;
imag_v_12 = 4 * t_12 * sa * ( ca - cb ) ...
    + 4 * u_12 * s2a * ( c2a - c2b );
v_12 = real_v_12 + 1i * imag_v_12;

v_22 = eps_2 + ( 3 * t_11 + t_22 ) * ca * cb + 2 * t_22 * c2a ...
    + 2 * r_11 * ( 2 * c3a * cb + c2b ) ...
    + 2 / sqrt(3) * r_12 * ( 4 * c3a * cb - c2b ) ...
    + ( 3 * u_11 + u_22 ) * c2a * c2b + 2 * u_22 * c4a;

H_TNN = [v_0, v_1, v_2; v_1', v_11, v_12; v_2', v_12', v_22];
