% 3-Band Tight-Binding MoS2
% Nach: Liu et el., Phys. Rev. B 88, 085433 (2013)
% Dritt-nächste Nachbarn in Symmetrie D_3h (d.h. in den Integralen ist das
% Hopping über die Schwefelatome durch die Symmetrie berücksichtigt).

function [ Ek, coeff, gradA_H, gradB_H ] = tb_MoS2_liu( params, kpts )

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
% E12     = params(21);
% Em10    = params(22);
% E10     = params(23);

k   = kpts;
sz  = size(k,2);
tri = size(k,3);

Ek      = zeros(6,sz,tri);
coeff   = complex(zeros(6,6,sz,tri));
gradA_H = complex(zeros(6,6,sz,tri));
gradB_H = complex(zeros(6,6,sz,tri));

%TB-Matrix
for tri_idx = 1:tri
    
    for kk = 1:sz
        
        H = complex(zeros(3,3));
        dA_H = complex(zeros(3,3));
        dB_H = complex(zeros(3,3));
        
        ktmp = [k(1,kk,tri_idx) k(2,kk,tri_idx)];
        A = ktmp(1) / 2;
        B = ktmp(2) / 2 * sqrt(3);
%         C = ktmp(3);
        
        %Diagonale
        %z2matlab
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
        % H_soc = kron(H,eye(2)) + lambda/2*kron(LS,[1 0;0 -1]);
        LS = [0 0 0;0 0 2*1i; 0 -2*1i 0];
        H1 = H + lambda/2.*LS;
        [V1, D1] = eig(H1*1e3);
        Ek1 = (real(diag(D1)));
        
        
        if kk == 1 && tri_idx == 1
            V1
            %             H2
        end
        
        
        H2 = H - lambda/2.*LS;
        [V2, D2] = eig(H2*1e3);
        Ek2 = (real(diag(D2)));
        
        Ek(1:3,kk,tri_idx) = Ek1(:);
        Ek(4:6,kk,tri_idx) = Ek2(:);
        coeff(1:3,1:3,kk,tri_idx) = V1;
        coeff(4:6,4:6,kk,tri_idx) = V2;
        

        % Ableitung
        %Diagonale
        %z2 ab=1
        dA_H(1,1) = -4*t0* ( sin(A)*cos(B) + sin(2*A) ) + ...
            -12*r0*( sin(3*A)*cos(B) )  ...
            -8*u0* ( sin(2*A)*cos(2*B) + sin(4*A) );
        
        dB_H(1,1) =  2*t0*( -2*cos(A)*sin(B) ) ...
            -4*r0*(  cos(3*A)*sin(B) + sin(2*B) ) + ...
            2*u0*( -4*cos(2*A)*sin(2*B) );
        
        
        %xy ab=4
        dA_H(2,2) = - ( t11 + 3*t22 )*sin(A)*cos(B) - 4*t11*sin( 2*A ) ...
            - 12*r11*sin( 3*A )*cos(B) ...
            - 2*( u11 + 3*u22 )*sin( 2*A )*cos(2*B) - 8*u11*sin( 4*A );
        
        dB_H(2,2) = - ( t11 + 3*t22 )*cos(A)*sin(B) ...
            - 4*r11*cos( 3*A )*sin(B) - 4*( r11 + sqrt(3)*r12 )*sin( 2*B ) ...
            - 2*( u11 + 3*u22 )*cos( 2*A )*sin( 2*B );
        
        
        %x2 ab=6
        dA_H(3,3) = - ( 3*t11 + t22 )*sin(A)*cos(B) - 4*t22*sin( 2*A ) + ...
            - 12*r11*sin( 3*A )*cos(B) ...
            - 8*sqrt(3)*r12*sin( 3*A )*cos(B) ...
            - 2*( 3*u11 + u22 )*sin( 2*A )*cos( 2*B ) - 8*u22*sin( 4*A );
        
        dB_H(3,3) = - ( 3*t11 + t22 )*cos(A)*sin(B) ...
            - 4*r11*( cos( 3*A )*sin(B)+sin( 2*B ) ) + ...
            4/sqrt(3)*r12*( sin( 2*B ) - 2*cos( 3*A )*sin(B) ) + ...
            - 2*( 3*u11 + u22)*cos( 2*A )*sin( 2*B );
        
        % Nebendiagonalterme
        % ab=2
        dA_H(1,2) = - 2*sqrt(3)*t2*cos(A)*sin(B) + 6*( r1 + r2 )*cos( 3*A )*sin(B) ...
            - 4*sqrt(3)*u2*cos( 2*A )*sin( 2*B ) ...
            +  1i*( 2*t1*cos(A)*( 2*cos(A) + cos(B) ) - 4*t1*sin(A).^2 + 6*( r1 - r2 )*cos( 3*A )*cos(B) ...
            + 4*u1*cos( 2.0*A )*( 2*cos( 2*A ) + cos( 2*B ) ) - 8*u1*sin( 2*A ).^2 );
        
        dB_H(1,2) = - 2*sqrt(3)*t2*sin(A)*cos(B)+2*(r1+r2)*sin(3*A)*cos(B) ...
            - 4*sqrt(3)*u2*sin(2*A)*cos(2*B) ...
            +  1i*( - 2*t1*sin(A)*sin(B) - 2*( r1 - r2 )*sin( 3*A )*sin(B) ...
            - 4*u1*sin( 2*A )*sin( 2*B ) );
        
        
        % ab=3
        dA_H(1,3) = 2*t2*( - 2*sin( 2*A ) + sin(A)*cos(B) ) ...
            - 2/sqrt(3)*( r1 + r2 )*( - 3*sin( 3*A )*cos(B) ) ...
            + 2*u2*( - 4*sin( 4*A ) + 2*sin( 2*A )*cos( 2*B ) ) ...
            +  1i*( - 2*sqrt(3)*t1*sin(A)*sin(B) ...
            + 2/sqrt(3)*sin(B)*( r1 - r2 )*( -3*sin( 3*A ) ) ...
            - 4*sqrt(3)*u1*sin( 2*A )*sin( 2*B ) );
        
        dB_H(1,3) = 2*t2*( cos(A)*sin(B) ) ...
            - 2/sqrt(3)*( r1 + r2 )*( - cos( 3*A )*sin(B) + 2*sin( 2*B ) ) ...
            + 4*u2*cos( 2*A )*sin( 2*B ) ...
            +  1i*( 2*sqrt(3)*t1*cos(A)*cos(B) ...
            + 2/sqrt(3)*( r1 - r2 )*( cos(B)*cos( 3*A ) + 2*cos(B).^2 - 2*sin(B).^2) ...
            + 4*sqrt(3)*u1*cos( 2*A )*cos( 2*B ) );
        
        
        % ab=5
        dA_H(2,3) =   sqrt(3)*( t22 - t11 )*cos(A)*sin(B) + 12*r12*cos( 3*A )*sin(B) ...
            + 2*sqrt(3)*( u22 - u11 )*cos( 2*A )*sin( 2*B ) ...
            +  1i*( 4*t12*cos(A)*( cos(A) - cos(B) ) - 4.0*t12*sin(A).^2  ...
            + 8*u12*cos( 2*A )*( cos( 2*A ) - cos( 2*B ) ) -8*u12*sin(2*A).^2);
        
        dB_H(2,3) =   sqrt(3)*( t22 - t11 )*sin(A)*cos(B) + 4*r12*sin( 3*A )*cos(B) ...
            + 2*sqrt(3)*( u22 - u11 )*sin( 2*A )*cos( 2*B ) ...
            +  1i*( 4*t12*sin(A)*sin(B) ...
            + 8*u12*sin( 2*A )*sin( 2*B ) );
        
        % ab=7
        dA_H(2,1) = conj( dA_H(1,2) );
        dB_H(2,1) = conj( dB_H(1,2) );
        
        % ab=8
        dA_H(3,1) = conj( dA_H(1,3) );
        dB_H(3,1) = conj( dB_H(1,3) );
        
        % ab=9
        dA_H(3,2) = conj( dA_H(2,3) );
        dB_H(3,2) = conj( dB_H(2,3) );
        
        gradA_H(1:3,1:3,kk,tri_idx) = dA_H; % * a / 2
        gradA_H(4:6,4:6,kk,tri_idx) = dA_H; % * a / 2
        
        gradB_H(1:3,1:3,kk,tri_idx) = dB_H; % * sqrt(3) / 2 * a
        gradB_H(4:6,4:6,kk,tri_idx) = dB_H; % * sqrt(3) / 2 * a
        
    end
    
end

Ek(1:3,:,:) = Ek(1:3,:,:) - max( max( Ek(1,:,:) ) );
Ek(4:6,:,:) = Ek(4:6,:,:) - max( max( Ek(1,:,:) ) );

return;

end