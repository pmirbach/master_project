function [pts] = pts_triangle_fun(pts, corners, epsn)
% Reihenfolge der Ecken: Unten links, Unten rechts, Oben

% Berechnung der Geradengleichungen
m = zeros(1,3);
b = m;
m(1) = (corners(2,2) - corners(2,1)) / (corners(1,2) - corners(1,1));
m(2) = (corners(2,3) - corners(2,1)) / (corners(1,3) - corners(1,1));
m(3) = (corners(2,3) - corners(2,2)) / (corners(1,3) - corners(1,2));
b(1) = corners(2,1) - m(1) * corners(1,1);
b(2) = corners(2,1) - m(2) * corners(1,1);
b(3) = corners(2,2) - m(3) * corners(1,2);

% Reduktion auf Punkte in Dreieck
a = pts(1,:) < min(corners(1,:)) - epsn ...
    | pts(1,:) > max(corners(1,:)) + epsn;
a = a | pts(2,:) < b(1) + m(1) * pts(1,:) - epsn;
a = a | pts(2,:) > b(2) + m(2) * pts(1,:) + epsn;
a = a | pts(2,:) > b(3) + m(3) * pts(1,:) + epsn;
pts(:,a) = [];

% Bestimmung der Gewichte
pts = [pts; 6 * ones(1,size(pts,2))];
for n = 1:3
    a = pts(2,:) > b(n) + m(n) * pts(1,:) - epsn ...
        & pts(2,:) < b(n) + m(n) * pts(1,:) + epsn;
    pts(3,a) = floor(pts(3,a) / 2);
end