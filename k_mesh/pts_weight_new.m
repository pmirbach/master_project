function [w] = pts_weight_new(pts, corners, epsn)

% Reihenfolge der Ecken: K, Gamma, K'

m(1) = corners(2,1) / corners(1,1);        % Steigung von Gamma zu K
m(2) = - m(1);


counter = zeros(1, size(pts,2));

for n = 1:2
    a = pts(2,:) > m(n) * pts(1,:) - epsn & pts(2,:) < m(n) * pts(1,:) + epsn;
    counter(a) = counter(a) + 1;
end

a = ( pts(1,:) > corners(1,1) - epsn & pts(1,:) < corners(1,1) + epsn );
counter(a) = counter(a) + 1;


w = 6 * ones(1, size(pts,2));

w( counter == 1 ) = 3;
w( counter == 2 ) = 1;