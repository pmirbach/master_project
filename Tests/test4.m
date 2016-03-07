clear variables
close all
clc

x = ones(1,10);
y = 1:6;
[X,Y] = meshgrid(x,y);

Y

W = diag([-1 1 1 -1 1 1])

W*Y