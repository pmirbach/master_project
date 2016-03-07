clear variables
close all
clc



fig1.h = test2;

ax = findobj(fig1.h(:),'type','axes');
hl = findobj(fig1.h(1),'type','legend');

ax(1).XLabel.String = 'ax1';
