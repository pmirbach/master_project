clear variables
close all
clc

% x = 0:0.01:10;
% y = -(x-5).^2 + 2;
% % plot(x,y)
% a = y < 0;
% 
% y(a) = nan;
% 
% % figure
% % plot(x,y)
% % 
% % axis([2 8 -0.5 2.5])
% 
% b = area(x,y);

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'normalized';
fig.Position = [0.1 0.1 0.8 0.8];

ax = gobjects(1,6);
for ii = 1:6
    ax(ii) = subplot(2,3,ii);
%     axis off
    colorbar
    ax(ii).Position = ax(ii).Position +[-.03 -.03 0 0];
end










