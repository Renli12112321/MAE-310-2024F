clear all; clc;

E      = 1e9;   %杨氏模量
nu     = 0.3;   %泊松比
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
mu     = E / (2 * (1 + nu));
L = 4;
Tx = 1e4;
R = 0.5;
r = @(x,y) sqrt((x+L/2)^2 + (y+L/2)^2);


% 定义极坐标下的变量函数
sigma_rr = @(x, y) Tx / 2 * (1 - R^2 / r(x, y)^2) + ...
    Tx / 2 * (1 - 4 * R^2 / r(x, y)^2 + 3 * R^4 / r(x, y)^4) * cos(2 * atan2(y + L / 2, x + L / 2));

sigma_tt = @(x, y) Tx / 2 * (1 + R^2 / r(x, y)^2) - ...
    Tx / 2 * (1 + 3 * R^4 / r(x, y)^4) * cos(2 * atan2(y + L / 2, x + L / 2));

sigma_rt = @(x, y) -Tx / 2 * (1 + 2 * R^2 / r(x, y)^2 - 3 * R^4 / r(x, y)^4) * sin(2 * atan2(y + L / 2, x + L / 2));

% 将极坐标变量转换为直角坐标变量
sigma_xy = @(x, y) sin(-atan2(y + L / 2, x + L / 2)) * cos(-atan2(y + L / 2, x + L / 2)) * (sigma_tt(x, y) - sigma_rr(x, y)) + ...
    sigma_rt(x, y) * (cos(-atan2(y + L / 2, x + L / 2))^2 - sin(-atan2(y + L / 2, x + L / 2))^2);

sigma_yy = @(x, y) sin(-atan2(y + L / 2, x + L / 2))^2 * sigma_rr(x, y) + ...
    cos(-atan2(y + L / 2, x + L / 2))^2 * sigma_tt(x, y) - ...
    2 * sin(-atan2(y + L / 2, x + L / 2)) * cos(-atan2(y + L / 2, x + L / 2)) * sigma_rt(x, y);

sigma_xx = @(x, y) cos(-atan2(y + L / 2, x + L / 2))^2 * sigma_rr(x, y) + ...
    sin(-atan2(y + L / 2, x + L / 2))^2 * sigma_tt(x, y) + ...
    2 * sin(-atan2(y + L / 2, x + L / 2)) * cos(-atan2(y + L / 2, x + L / 2)) * sigma_rt(x, y);
%D      = [lambda + 2 * mu, lambda,          0;
%lambda,          lambda + 2 * mu, 0;
%0,               0,               mu];%根据note L11算出D矩阵

%不是，为什么L11开头的公式直接就是错的，两种方法算出来的D不相等啊

D      = [E/(1-nu^2),        nu*E/(1-nu^2),          0;
    nu*E/(1-nu^2),    E/(1-nu^2),              0;
    0,                0,               E/2/(1+nu)];

% 体力分量
f_x = @(x, y) ...
    (E * nu * ((x - 1) * (y - 1) + x * y + x * (y - 1) + y * (x - 1))) / (nu^2 - 1) ...
    + (2 * E * y * (y - 1) - E * (nu / 2 - 1 / 2) * (x - 1) * (y - 1) - E * (nu / 2 - 1 / 2) * (x * y + 2 * x * (x - 1) + x * (y - 1) + y * (x - 1))) / (nu^2 - 1);

f_y = @(x, y) ...
    (E * nu * ((x - 1) * (y - 1) + x * y + x * (y - 1) + y * (x - 1))) / (nu^2 - 1) ...
    + (2 * E * x * (x - 1) - E * (nu / 2 - 1 / 2) * (x - 1) * (y - 1) - E * (nu / 2 - 1 / 2) * (x * y + x * (y - 1) + y * (x - 1) + 2 * y * (y - 1))) / (nu^2 - 1);

% 体力矢量
body_force = @(x, y) [f_x(x, y); f_y(x, y)];

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


% mesh generation
run('quarter_plate_with_hole_quad2.m')
[x_coor, y_coor] = deal(msh.POS(:, 1), msh.POS(:, 2));  % 解构坐标列
n_el = size(msh.QUADS, 1);  % total number of elements
n_np = msh.nbNod;  % total number of nodal points
line_in = msh.LINES;  % 线信息
n_en   = 4;               % number of nodes in an element


% IEN array
%IEN = zeros(n_el, n_en);
IEN = msh.QUADS(:,1:4);

%for ex = 1 : n_el_x
    %for ey = 1 : n_el_y
        %ee = (ey-1) * n_el_x + ex; % element index
        %IEN(ee, 1) = (ey-1) * n_np_x + ex;
        %IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
        %IEN(ee, 3) =  ey    * n_np_x + ex + 1;
        %IEN(ee, 4) =  ey    * n_np_x + ex;
    %end
%end


% ID array
ID = zeros(n_np,2);
counter = 0;
for ny = 2 : n_np_y - 1
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;
        counter = counter + 1;
        ID(index,1) = counter;
        counter = counter + 1;
        ID(index,2) = counter;
    end
end

n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );

    k_ele = zeros(2*n_en, 2*n_en); % element stiffness matrix
    f_ele = zeros(2*n_en, 1);    % element load vector

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        bf = body_force(x_l, y_l);
        Ba = zeros(3, 2); % 应变-位移矩阵

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            for dir = 1:2 % 方向：1=u, 2=v
                f_ele(2*aa-2+dir) = f_ele(2*aa-2+dir) + weight(ll) * detJ * bf(dir) * Na;
            end

            Ba(1, 1) = Na_x;
            Ba(2, 2) = Na_y;
            Ba(3, 1) = Na_y;
            Ba(3, 2) = Na_x;

            Bb = zeros(3,2); % 应变-位移矩阵
            for bb = 1 : n_en
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                Bb(1, 1) = Nb_x;
                Bb(2, 2) = Nb_y;
                Bb(3, 1) = Nb_y;
                Bb(3, 2) = Nb_x;

                for ii = 1:2
                    for jj = 1:2
                        k_ele(2*(aa-1)+ii,2*(bb-1)+jj) = k_ele(2*(aa-1)+ii,2*(bb-1)+jj)+...
                            weight(ll)*detJ*((ii==1)*[1,0]+(ii==2)*[0,1])* Ba'*D* Bb*...
                            ((jj==1)*[1;0]+(jj==2)*[0;1]);
                    end
                end

                %k_ele = k_ele + Ba' * D * Bb * detJ * weight(ll);

            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for aa = 1:n_en
        for dirA = 1:2 % 方向：1=u, 2=v
            PP = ID (IEN(ee , aa), dirA);
            if PP > 0
                F(PP) = F(PP) + f_ele(2*aa-2+dirA);
                for bb = 1:n_en
                    for dirB = 1:2 % 方向：1=u, 2=v
                        QQ = ID (IEN(ee , bb), dirB);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(2*aa-2+dirA, 2*bb-2+dirB);
                        end
                    end
                end
            end
        end
    end


end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 2);



% EOF