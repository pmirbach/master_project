function [x] = red_to_BZ(x)

% Definition der Symmetrieoperationen
C3 = [cos(2 * pi / 3), -sin(2 * pi / 3), 0; ...
    sin(2 * pi / 3) cos(2 * pi / 3), 0; 0, 0, 1];
sigma1 = [1, 0, 0 ; 0 ,-1,0 ; 0, 0, 1];

% Bei 2 dim Vektoren - Reduktion der Matrizen
if size(x,1) < 3
    C3 = C3(1:2,1:2);
    sigma1 = sigma1(1:2,1:2);
end

% Anwenden der Symmetrieoperationen
x(:,:,3) = (C3 * x(:,:,1));
x(:,:,5) = (C3 * x(:,:,3));
x(:,:,6) = (sigma1 * x(:,:,1));
x(:,:,2) = (C3 * x(:,:,6));
x(:,:,4) = (C3 * x(:,:,2));