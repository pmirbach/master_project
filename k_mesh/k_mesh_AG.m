function k_mesh_AG(Ctrl, Para)

Ctrl.plot.k_meshAG = [ 1 1 1 1 ];

a = 0.318;

b1 = 2 * pi / a * [ 2 / sqrt(3) ; 0 ];
b2 = 2 * pi / a * [ 1 / sqrt(3) ; 1 ];

newsymmpts = 2 * pi / ( 3 * a ) * [ [ sqrt(3); 1 ] , [ 0; 0 ], [ sqrt(3); -1 ] ];

% Density of k-mesh
qr = Ctrl.k_mesh_mp.qr;

% Ascending and descending matrices for rhomboid k-mesh including red. BZ
ur = ( -2 * qr / 3 : 2 * qr / 3 ) / qr;
A_asc = ones( numel(ur),1 ) * ur;
A_desc = -ones(numel(ur),1) * ur;

% Ascending and descending x component of k-mesh
k_rhomboid_x_asc = A_asc' * b1(1) + A_asc * b2(1);
k_rhomboid_x_desc = A_desc' * b1(1) + A_desc * b2(1);
% y component of k-mesh
k_rhomboid_y = A_asc' * b1(2) + A_asc * b2(2);

% All k points (asc. and desc.)
k_all_asc = [k_rhomboid_x_asc(:), k_rhomboid_y(:)]';
k_all_desc = [k_rhomboid_x_desc(:), k_rhomboid_y(:)]';

% Classification of k points above, below and on bisectrix in red. BZ
k_rhomboid_below = k_all_asc( : , k_all_asc(2,:) < 0 );
k_rhomboid_bisec = k_all_asc( : , k_all_asc(2,:) == 0 );
k_rhomboid_above = k_all_desc( : , k_all_desc(2,:) > 0 );

% Flip k points to adress high-symmetrie points with index 1 and end
k_rhomboid_below = k_rhomboid_below(: , end:-1:1);
k_rhomboid_above = k_rhomboid_above(: , end:-1:1);

% All k points in rhomboid mesh with desired indexing
k_rhomboid = [ k_rhomboid_above , k_rhomboid_bisec , k_rhomboid_below ];

% Exclude all k points outside red. BZ
[ k_redBZ , wk ] = kpts_in_redBZ( k_rhomboid , newsymmpts );

% Find equivalent indexes of upper, lower and middle k points
[ Nrk , ind ] = identify_indexes( k_redBZ );

% Create k points in entire BZ
[ k_BZ ] = create_k_in_BZ( k_redBZ );


% Create k-mesh test plots if wanted
% k_mesh_test_plots( Ctrl , k_redBZ , wk , ind , k_BZ )

end


function [ pts , weight ] = kpts_in_redBZ( pts , corners )

% Order of corners: K, Gamma, K'

% Numeric margin of error
epsn = 30 * eps;

% Determine the slope for lines from Gamma to K / K'
m = corners(2,1) / corners(1,1);

% Find kpts outside rectangle around red. BZ
bool_exclude = pts(1,:) < min(corners(1,:)) - epsn | pts(1,:) > max(corners(1,:)) + epsn;
bool_exclude = bool_exclude | pts(2,:) < min(corners(2,:)) - epsn | pts(2,:) > max(corners(2,:)) + epsn;
% Find kpts above line Gamma to K and below line Gamma to K'
bool_exclude = bool_exclude | pts(2,:) > m * pts(1,:) + epsn | pts(2,:) < -m * pts(1,:) - epsn;
% Exclude the found k points
pts(:,bool_exclude) = [];


% Determination of the integer weight frequency in entire BZ
% High-symmetrie points:    1
% Pts on edge of red. BZ:   3
% Remaining points:         6

m = [ m , -m ];
% Counter for points on edges - on 2 edges -> corner point
counter = zeros(1, size(pts,2));
% Find points on one of the three edges
for ni = 1:2
    bool_on_edge = pts(2,:) > m(ni) * pts(1,:) - epsn & pts(2,:) < m(ni) * pts(1,:) + epsn;
    counter(bool_on_edge) = counter(bool_on_edge) + 1;
end
bool_on_edge = ( pts(1,:) > corners(1,1) - epsn & pts(1,:) < corners(1,1) + epsn );
counter(bool_on_edge) = counter(bool_on_edge) + 1;

weight = 6 * ones(1, size(pts,2));

weight( counter == 1 ) = 3;
weight( counter == 2 ) = 1;

end


function [ Nrk , ind ] = identify_indexes( k_redBZ )

Nrk = size(k_redBZ,2);
Nrmitte = sum( k_redBZ(2,:) == 0 );
NrSide = ( Nrk - Nrmitte ) / 2;

ind.up = 1 : NrSide;
ind.mid = NrSide + 1 : NrSide + Nrmitte;
ind.dwn = Nrk : -1 : Nrk - NrSide + 1;

end


function [ k_BZ ] = create_k_in_BZ( k_redBZ )

k_BZ = zeros( [size(k_redBZ),6] );
k_BZ(:,:,1) = k_redBZ;

% Counterpart of tri 1 is tri 4: Same ky, mirrored kx
k_BZ(:,:,4) = k_BZ(:,:,1);
k_BZ(1,:,4) = - k_BZ(1,:,4);

% Other triangles through Rotation by 60 degrees
C3 = [ cos(2 * pi / 3), -sin(2 * pi / 3); ...
    sin(2 * pi / 3) cos(2 * pi / 3) ];

k_BZ(:,:,3) = C3 * k_BZ(:,:,1);
k_BZ(:,:,5) = C3 * k_BZ(:,:,3);
k_BZ(:,:,6) = C3 * k_BZ(:,:,4);
k_BZ(:,:,2) = C3 * k_BZ(:,:,6);

end


function k_mesh_test_plots( Ctrl , k_redBZ , wk , ind , k_BZ )


if Ctrl.plot.k_meshAG(1) == 1
    a
end

if Ctrl.plot.k_meshAG(2) == 1
    a
end

if Ctrl.plot.k_meshAG(3) == 1
    a
end

if Ctrl.plot.k_meshAG(4) == 1
    a
end


ind = 30;
figure; hold on
for ii = 1:6
    plot(k_BZ(1,:,ii),k_BZ(2,:,ii),'x')
    plot(k_BZ(1,ind,ii),k_BZ(2,ind,ii),'ko')
end



% Test Plot: All in One

kind = 1:Nrk;

bool_or = false(3, Nrk);
bool_or(1,:) = kind <= max( ind_oben );
bool_or(2,:) = kind > max( ind_oben ) & kind < min( ind_unten );
bool_or(3,:) = kind >= min( ind_unten );

bool_w = false(3, Nrk);
bool_w(1,:) = wk == 6;
bool_w(2,:) = wk == 3;
bool_w(3,:) = wk == 1;

colors = {'b','r','g','c','m'};
marker = {'h','x','^','+','v'};
figure; hold on
for ni = 1:3
    for nj = 1:3
        asdf = [ colors{nj}, marker{ni} ];
        plot(k(1,bool_w(ni,:) & bool_or(nj,:),1),k(2,bool_w(ni,:) & bool_or(nj,:),1),asdf)
    end
end



% Test Plot: Gewichtung
in = wk == 6;
on = wk == 3;
symm = wk == 1;
figure; hold on
plot(kall(1,in,2),kall(2,in,2),'rx')
plot(kall(1,on,2),kall(2,on,2),'r^')
plot(kall(1,symm,2),kall(2,symm,2),'rs')

% Test Plot: Indizierung 1.1 - Speicherordnung
figure; hold on
for ii = 1:6
    plot(kall(1,:,ii),kall(2,:,ii),'x')
end
% Test Plot: Indizierung 1.2 - aequivalente Punkte
hold on
ind = 1;
ind_eq = Nrk - ind + 1;
plot(kall(1,ind,1),kall(2,ind,1),'ko')
plot(kall(1,ind_eq,1),kall(2,ind_eq,1),'bo')


% Test Plot: Indizierung 3 - Indizes unten, mitte, oben
color_matrix = [ 160 82 45 ; 205 133 63 ; 222 184 135 ; 85 107 47 ; 107 142 35 ; 189 183 107] / 255;
opac = 1/10;

figure; hold on
for ni = 1:6
    scatter(kall(1,ind_oben,ni),kall(2,ind_oben,ni),'^','MarkerEdgeColor',color_matrix(ni,:),...
        'MarkerFaceColor',color_matrix(ni,:),'MarkerFaceAlpha',opac)
    scatter(kall(1,ind_mitte,ni),kall(2,ind_mitte,ni),'s','MarkerEdgeColor',color_matrix(ni,:),...
        'MarkerFaceColor',color_matrix(ni,:),'MarkerFaceAlpha',opac)
    scatter(kall(1,ind_unten,ni),kall(2,ind_unten,ni),'v','MarkerEdgeColor',color_matrix(ni,:),...
        'MarkerFaceColor',color_matrix(ni,:),'MarkerFaceAlpha',opac)
end

testinds = [121, 250];

for ii = 1:length(testinds)
    
    for ni = 1:6
        
        indup = ind_oben( testinds( ii ) );
        inddwn = ind_unten( testinds( ii ) );
        
        plot(kall(1,indup,ni),kall(2,indup,ni),'ko')
        plot(kall(1,inddwn,ni),kall(2,inddwn,ni),'ro')
        
    end
    
end








% Test Plot: Entire BZ
figure; hold on
for ii = 1:6
    plot(kall(1,:,ii),kall(2,:,ii),'x')
end


end
