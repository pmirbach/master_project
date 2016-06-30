clear variables
% close all
clc

wahl = 4;

a = 0.319;

switch wahl
    case 1
        a1 = a / 2 * [ 1 , sqrt(3) , 0 ];
        a2 = a / 2 * [ 1 , -sqrt(3) , 0 ];
    case 2
        a1 = a / 2 * [ 1 , sqrt(3) , 0 ];
        a2 = a / 2 * [ -1 , sqrt(3) , 0 ];
    case 3
        a1 = a / 2 * [ 2 , 0 , 0 ];
        a2 = a / 2 * [ -1 , sqrt(3) , 0 ];
    case 4
        a1 = a / 2 * [ sqrt(3) , 1 , 0 ];
        a2 = a / 2 * [ 0 , 2 , 0 ];
    case 5
        a1 = a / 2 * [ sqrt(3) , 1 , 0 ];
        a2 = a / 2 * [ 0 , 2 , 0 ];
end


a3 = [ 0 , 0 , 1 ];

V = dot( a1 , cross( a2 , a3 ) );

b1 = 2 * pi * cross( a2 , a3 ) / V;
b2 = 2 * pi * cross( a3 , a1 ) / V;


figure
subplot(1,2,1)
quiver(0,0,a1(1),a1(2))
hold on
quiver(0,0,a2(1),a2(2))
dd = max( abs( a1 ) ) + 0.1;
axis([-dd dd -dd dd])
axis equal


subplot(1,2,2)
quiver(0,0,b1(1),b1(2))
hold on
quiver(0,0,b2(1),b2(2))
ee = max( abs( b1 ) ) + 0.1;
% axis([-ee ee -ee ee])


% quiver(0,0,-b1(1),-b1(2))
% quiver(0,0,-b2(1),-b2(2))
% 
% if norm( b1 + b2 ) < norm( b1 - b2 )
%     quiver(0,0,b1(1)+b2(1),b1(2)+b2(2))
%     quiver(0,0,-b1(1)-b2(1),-b1(2)-b2(2))
% else
%     quiver(0,0,b1(1)-b2(1),b1(2)-b2(2))
%     quiver(0,0,-b1(1)+b2(1),-b1(2)+b2(2))    
% end



axis equal

%%
% % Generate hexagonal grid
% Rad3Over2 = sqrt(3) / 2;
% [X Y] = meshgrid(0:1:41);
% n = size(X,1);
% X = Rad3Over2 * X;
% Y = Y + repmat([0 0.5],[n,n/2]);
% 
% % Plot the hexagonal mesh, including cell borders
% [XV YV] = voronoi(X(:),Y(:)); plot(XV,YV,'k-')
% axis equal, axis([10 20 10 20]), zoom on










