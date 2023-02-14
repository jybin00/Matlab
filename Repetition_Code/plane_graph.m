clear
close all
clc

A = 2;
B = 2;
C = 3;
D = 3;

[x y] = meshgrid(-1.5:0.1:1.5, -1.5:0.1:1.5);

z = -x -y;

surf(x, y, z)
hold on
plot3(1, 1, 1, "o")
hold on
plot3(-1, -1, -1, "o")

fun = @(x,y,z) 1/(2*pi)^(3/2)*exp(-((x).^2 + (y).^2 + (z).^2)/2);

xmin = -sqrt(3);
xmax = sqrt(3);
ymin = @(x) -sqrt(3-x.^2);
ymax = @(x) sqrt(3-x.^2);
zmin = @(x,y) -sqrt(3-x.^2 - y.^2);
zmax = @(x,y) sqrt(3-x.^2 - y.^2);

q = integral3(fun, xmin, xmax, ymin, ymax, zmin, zmax, 'Method', 'tiled')