clear all; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x)
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 2;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
for n_el = 2:2:16      % mesh with element number from 2 to 16
    hh = 1 / n_el /n_en;
    x_coor = 0 : hh : 1;
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations

%Setup the IEN
IEN = zeros(n_el, n_en);
for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1) * pp + aa;
  end
end

% Setup the quadrature rule
n_int = 10;
[xi, weight] = Gauss(n_int, -1, 1);
