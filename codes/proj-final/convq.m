clear all; clc;

E      = 1e9;   %杨氏模量
nu     = 0.3;   %泊松比
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
mu     = E / (2 * (1 + nu));
D      = [lambda + 2 * mu, lambda,          0;
    lambda,          lambda + 2 * mu, 0;
    0,               0,               mu];%根据note L11算出D矩阵


%loop修改思路:
%应力应变问题每个节点有两个自由度（位移u,v）
%应力应变问题中，刚度矩阵由应变-位移矩阵B和材料刚度矩阵D乘积确定
%形函数导数用于构建应变-位移矩阵B，而非计算梯度
%设定精确的位移，一阶导确定应变精确解，结合D矩阵可得应力精确解
%

% 人工设定位移解
u_exact = @(x, y) x*(1-x)*y*(1-y);
v_exact = @(x, y) x*(1-x)*y*(1-y);

% 一阶导数（得出应变的精确解）
u_x = @(x, y) (1 - 2 * x) * y * (1 - y); % du/dx
u_y = @(x, y) x * (1 - x) * (1 - 2 * y); % du/dy
v_x = @(x, y) (1 - 2 * x) * y * (1 - y); % dv/dx
v_y = @(x, y) x * (1 - x) * (1 - 2 * y); % dv/dy

%应变精确解
epsilon_xx = @(x, y) u_x(x, y);
epsilon_yy = @(x, y) v_y(x, y);
gamma_xy = @(x, y) u_y(x, y) + v_x(x, y);

% 应力精确解
sigma_xx = @(x, y) (lambda + 2 * mu) * epsilon_xx(x, y) + lambda * epsilon_yy(x, y);
sigma_yy = @(x, y) (lambda + 2 * mu) * epsilon_yy(x, y) + lambda * epsilon_xx(x, y);
sigma_xy = @(x, y) mu * gamma_xy(x, y);

% 二阶导数
u_xx = @(x, y) -2 * y * (1 - y); % d²u/dx²
u_yy = @(x, y) -2 * x * (1 - x); % d²u/dy²
v_xx = @(x, y) -2 * y * (1 - y); % d²v/dx²
v_yy = @(x, y) -2 * x * (1 - x); % d²v/dy²


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

% mesh size loop
logeL2 = zeros(8,1);
logeH1 = zeros(8,1);
logh   = zeros(8,1);
for n_me = 6:2:20

    % mesh generation
    n_en   = 4;               % number of nodes in an element
    n_el_x = n_me;               % number of elements in x-dir
    n_el_y = n_me;               % number of elements in y-dir
    n_el   = n_el_x * n_el_y; % total number of elements

    n_np_x = n_el_x + 1;      % number of nodal points in x-dir
    n_np_y = n_el_y + 1;      % number of nodal points in y-dir
    n_np   = n_np_x * n_np_y; % total number of nodal points

    x_coor = zeros(n_np, 1);
    y_coor = x_coor;

    hx = 1.0 / n_el_x;        % mesh size in x-dir
    hy = 1.0 / n_el_y;        % mesh size in y-dir

    % generate the nodal coordinates
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx; % nodal index
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end

    % IEN array
    IEN = zeros(n_el, n_en);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex; % element index
            IEN(ee, 1) = (ey-1) * n_np_x + ex;
            IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
            IEN(ee, 3) =  ey    * n_np_x + ex + 1;
            IEN(ee, 4) =  ey    * n_np_x + ex;
        end
    end


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

    for ii = 1 : n_np
        for jj = 1 : 2
            index = ID(ii , jj);
            if index > 0
                disp(ii , jj) = dn(index);
            else
                % modify disp with the g data. Here it does nothing because g is zero
            end
        end
    end
    save("disp", "disp", "n_el_x", "n_el_y");
    % calculate the error
    L2 = 0; H1 = 0;

    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        y_ele = y_coor(IEN(ee,:));
        u_ele = disp( IEN(ee, :), :);

        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            uh = zeros(2,1);
            uh_x = zeros(2,1); uh_y = zeros(2,1);
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

            for aa = 1 : n_en
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                for dir = 1:2
                    uh(dir)   = uh(dir)   + u_ele(aa,dir) * Na;
                    uh_x(dir) = uh_x(dir) + u_ele(aa,dir) * Na_x;
                    uh_y(dir) = uh_y(dir) + u_ele(aa,dir) * Na_y;
                end
            end

            L2 = L2 + weight(ll) * (uh - u_exact(x_l,y_l))'*(uh - u_exact(x_l,y_l)) * detJ;
            H1 = H1 + weight(ll) * ((uh_x - u_x(x_l,y_l) )'*(uh_x - u_x(x_l,y_l) ) + (uh_y - u_y(x_l,y_l) )'*(uh_y - u_y(x_l,y_l) )) * detJ;

        end
    end

    L2 = sqrt(L2);
    H1 = sqrt(H1);

    % 保存该mesh结果到数组
    logeL2(n_me/2-2) = log(L2);
    logeH1(n_me/2-2) = log(H1);
    logh  (n_me/2-2) = log(hx);

end
% save the solution vector and number of elements to disp with name
% HEAT.mat
%save("HEAT", "disp", "n_el_x", "n_el_y");
% 画图
plot(logh, logeL2, '-r','LineWidth',3);
xlabel('log(h)');
ylabel('log(error L2)');
hold on;
p1 = polyfit(logh,logeL2,1);
y1 = polyval(p1,logh);
plot(logh,y1, '-b','LineWidth',1);

figure
plot(logh, logeH1, '-k','LineWidth',3);
xlabel('log(h)');
ylabel('log(error H1)');
hold on;
p2 = polyfit(logh,logeH1,1);
y2 = polyval(p2,logh);
plot(logh,y2, '-b','LineWidth',1);

% EOF