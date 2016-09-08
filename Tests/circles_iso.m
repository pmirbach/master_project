clear variables
close all
clc





figure; hold on


radii = 1:7 ;
nr_r = numel(radii);

% Kreisgerüste
for ii = 1:nr_r
    rectangle('Position',[ -radii(ii) , -radii(ii) , 2 * radii(ii) , 2 * radii(ii) ],'Curvature',[ 1 , 1 ])
end
scatter(0,0,10,'k','MarkerFaceColor' , 'k')
axis equal



% Anzahl Kreise je Radius

nr_pts_iso_r = 1 : 6 : ( nr_r * 6 + 1 );

radius_pts = 0.5;

for ii = 1:nr_r
    
    rectangle('Position',[ -radii(ii) , -radii(ii) , 2 * radii(ii) , 2 * radii(ii) ],'Curvature',[ 1 , 1 ])
    
end



