function [w] = pts_weight(pts, corners, epsn)

% Berechnung der Geradengleichungen
m = zeros(1,3);
b = m;
m(1) = (corners(2,2) - corners(2,1)) / (corners(1,2) - corners(1,1));
m(2) = (corners(2,3) - corners(2,1)) / (corners(1,3) - corners(1,1));
m(3) = (corners(2,3) - corners(2,2)) / (corners(1,3) - corners(1,2));
b(1) = corners(2,1) - m(1) * corners(1,1);
b(2) = corners(2,1) - m(2) * corners(1,1);
b(3) = corners(2,2) - m(3) * corners(1,2);

w = 6 * ones(1, size(pts,2));
for n = 1:3
    a = pts(2,:) > b(n) + m(n) * pts(1,:) - epsn ...
        & pts(2,:) < b(n) + m(n) * pts(1,:) + epsn;
    w(a) = floor(w(a) / 2);
end