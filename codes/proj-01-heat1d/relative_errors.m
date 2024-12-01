clear all; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3;  % f(x)
g = 1.0;            % u    = g  at x = 1
h = 0.0;            % -u,x = h  at x = 0
ex_u = @(x) x.^5;   % exact solution
ex_u2= @(x) x.^10;  % exact u square
ex_du= @(x) 5*x.^4; % exact u,x
ex_d2= @(x) 25*x.^8;% exact u,x square
denoL2 = sqrt(integral(ex_u2, 0, 1)); % L2分母
denoH1 = sqrt(integral(ex_d2, 0, 1)); % H1分母

% Setup the quadrature rule
n_int = 10;
[xi, weight] = Gauss(n_int, -1, 1);

% 存储不同mesh结果的数组
logeL2 = zeros(8,1);
logeH1 = zeros(8,1);
logh   = zeros(8,1);

% Setup the mesh
pp   = 2;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
for n_el = 2:2:16      % mesh with element number from 2 to 16
    hh = 1 / n_el /n_en;
    x_coor = 0 : hh : 1;
    n_np = n_el * pp + 1;  % number of nodal points
    n_eq = n_np - 1;       % number of equations

    % Setup the IEN
    IEN = zeros(n_el, n_en);
    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end

    % Setup the ID array
    ID = 1 : n_np;
    ID(end) = 0;

    % Stiffness matrix
    K = zeros(n_eq, n_eq);
    F = zeros(n_eq, 1);

    % Assembly of the stiffness matrix and load vector
    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
        f_ele = zeros(n_en, 1);    % allocate a zero element load vector

        x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

        % quadrature loop
        for qua = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : n_en
                x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
            end
            dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : n_en
                f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
                end
            end
        end

        % Assembly of the matrix and vector based on the ID or LM data
        for aa = 1 : n_en
            P = ID(IEN(ee,aa));
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb));
                    if(Q > 0)
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                    end
                end
            end
        end
    end

    % ee = 1 F = NA(0)xh
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

    % Solve Kd = F equation
    d_temp = K \ F;
    disp = [d_temp; g];

    n_sam = 20;
    xi_sam = -1 : (2/n_sam) : 1;

    x_sam = zeros(n_el * n_sam + 1, 1);
    y_sam = x_sam; % store the exact solution value at sampling points
    u_sam = x_sam; % store the numerical solution value at sampling pts
    el2_sam = x_sam; % 存储L2差值
    
    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );

        if ee == n_el
            n_sam_end = n_sam+1;
        else
            n_sam_end = n_sam;
        end

        for ll = 1 : n_sam_end
            x_l = 0.0;
            u_l = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
            end

            x_sam( (ee-1)*n_sam + ll ) = x_l;
            u_sam( (ee-1)*n_sam + ll ) = u_l;
            y_sam( (ee-1)*n_sam + ll ) = x_l^5;
            el2_sam( (ee-1)*n_sam + ll ) = u_l-x_l^5;

        end
    end

    % 计算出eL2的分子部分
    numeL2 = 0;
    numeH1 = 0;
    for ee = 1 :  n_el * n_sam
        numeL2 = numeL2 + (x_sam(ee+1)-x_sam(ee))*el2_sam(ee)^2;
        numeH1 = numeH1 + (x_sam(ee+1)-x_sam(ee))*((u_sam(ee+1)-u_sam(ee))/(x_sam(ee+1)-x_sam(ee))-ex_du(x_sam(ee)))^2;
    end

    % 保存该mesh结果到数组
    numeL2 = sqrt(numeL2);
    numeH1 = sqrt(numeH1);
    logeL2(n_el/2) = log(numeL2/denoL2);
    logeH1(n_el/2) = log(numeH1/denoH1);
    logh  (n_el/2) = log(hh);

end

% 画图

plot(logh, logeL2, '-r','LineWidth',3);
xlabel('log(h)');
ylabel('log(error L2)');

hold on;

figure
plot(logh, logeH1, '-k','LineWidth',3);
xlabel('log(h)');
ylabel('log(error H1)');

