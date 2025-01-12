clear all; clc; clf;

load("disp.mat");

hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;

[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
hold on
Z = reshape(disp(:,1), n_np_x, n_np_y)';
surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);

% EOF