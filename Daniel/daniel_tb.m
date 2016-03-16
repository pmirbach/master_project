function [H1,H2] = daniel_tb(k, params)

a = 0.319;

E1  = params(1);
E2  = params(2);
t0  = params(3);
t1  = params(4);
t2  = params(5);
t11 = params(6);
t12 = params(7);
t22 = params(8);
r0  = params(9);
r1  = params(10);
r2  = params(11);
r11 = params(12);
r12 = params(13);
u0  = params(14);
u1  = params(15);
u2  = params(16);
u11 = params(17);
u12 = params(18);
u22 = params(19);

lambda  = params(20);

% A = k(1) * pi / 2;
% B = k(2) * pi / 2 * sqrt(3);

A = k(1) / 2 * a;
B = sqrt(3) * k(2) / 2 * a;


H = complex(zeros(3,3));

%Diagonale
%z2
H(1,1) = E1 + 2*t0*(2*cos(A)*cos(B)+cos(2*A)) + ...
    2*r0*(2*cos(3*A)*cos(B)+cos(2*B)) + ...
    2*u0*(2*cos(2*A)*cos(2*B)+cos(4*A));
%xy
H(2,2) = E2 + (t11+3*t22)*cos(A)*cos(B)+2*t11*cos(2*A) + ...
    4*r11*cos(3*A)*cos(B)+2*(r11+sqrt(3)*r12)*cos(2*B) + ...
    (u11+3*u22)*cos(2*A)*cos(2*B)+2*u11*cos(4*A);
%x2
H(3,3) = E2 + (3*t11+t22)*cos(A)*cos(B)+2*t22*cos(2*A) + ...
    2*r11*(2*cos(3*A)*cos(B)+cos(2*B)) + ...
    2/sqrt(3)*r12*(4*cos(3*A)*cos(B)-cos(2*B)) + ...
    (3*u11+u22)*cos(2*A)*cos(2*B)+2*u22*cos(4*A);


%# Nebendiagonalterme dd
H(1,2) = -2*sqrt(3)*t2*sin(A)*sin(B)+2*(r1+r2)*sin(3*A)*sin(B) - ...
    2*sqrt(3)*u2*sin(2*A)*sin(2*B) + ...
    1i*(2*t1*sin(A)*(2*cos(A)+cos(B))+2*(r1-r2)*sin(3*A)*cos(B) + ...
    2*u1*sin(2*A)*(2*cos(2*A)+cos(2*B)));

H(1,3) = 2*t2*(cos(2*A)-cos(A)*cos(B)) - ...
    2/sqrt(3)*(r1+r2)*(cos(3*A)*cos(B)-cos(2*B)) + ...
    2*u2*(cos(4*A)-cos(2*A)*cos(2*B)) + ...
    1i*(2*sqrt(3)*t1*cos(A)*sin(B) + ...
    2/sqrt(3)*sin(B)*(r1-r2)*(cos(3*A)+2*cos(B)) + ...
    2*sqrt(3)*u1*cos(2*A)*sin(2*B));

H(2,3) = sqrt(3)*(t22-t11)*sin(A)*sin(B)+4*r12*sin(3*A)*sin(B) + ...
    sqrt(3)*(u22-u11)*sin(2*A)*sin(2*B) + ...
    1i*(4*t12*sin(A)*(cos(A)-cos(B)) + ...
    4*u12*sin(2*A)*(cos(2*A)-cos(2*B)));


H(2,1) = conj(H(1,2));
H(3,1) = conj(H(1,3));
H(3,2) = conj(H(2,3));


% SOC
% H_soc = kron(H,eye(2)) + lambda/2*kron(LS,[1 0;0 -1]);;
LS = [0 0 0;0 0 2*1i; 0 -2*1i 0];


H1 = H + lambda/2.*LS;

H2 = H - lambda/2.*LS;