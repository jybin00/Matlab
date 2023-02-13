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