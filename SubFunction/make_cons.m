%% 根据约束条件，制作出约束矩阵

function [x0, lb, ub, A, b, F_phi_cons] = make_cons(F_data, phi_data, para, cons_data, control, F_phi_cons)
knots = para.knots;     % 记录了多项式的节点个数
order = para.order;     % 记录了多项式的阶次，下区间和上区间
rho = para.rho;         % 记录了电荷密度函数
range = para.range;     % 记录了F和phi函数的取值范围

alpha_fcc_gsf = para.alpha_fcc_gsf;     % 保存着FCC晶格广义平面缺陷能的几何因子
alpha_fcc = para.alpha_fcc;             % 保存着FCC平衡态的几何因子
N_fcc = para.N_fcc;
N_fcc_111 = para.N_fcc_111;
N_fcc_110 = para.N_fcc_110;
N_fcc_100 = para.N_fcc_100;

% 以后需要将这部分代码转移入结构体中
drho = (rho(3:end, 2) - rho(1:end-2, 2)) ./ (rho(3:end, 1) - rho(1:end-2, 1));
drho_x = rho(2:end-1, 1);

% 用于保存约束量
a_fcc = cons_data.lattice(1);           % FCC,HCP和BCC的晶格常数
a_hcp = cons_data.lattice(2);
c_hcp = cons_data.lattice(3);
a_bcc = cons_data.lattice(4);
fcc_elastic = cons_data.fcc_elastic;    % FCC,HCP和BCC的弹性常数
% hcp_elastic = cons_data.hcp_elastic;
% bcc_elastic = cons_data.bcc_elastic;
sub_energy = cons_data.sub_energy;      % FCC,HCP和BCC的升华能
sub_fcc = sub_energy(1);
sub_hcp = sub_energy(2);
sub_bcc = sub_energy(3);

% 用于保存约束条件
F_phi_range = control.F_phi_range;
fcc_relation_C12_C44 = control.fcc_relation_C12_C44;
dF_restrict = control.dF_restrict;
fcc_dF_restrict_equal = control.fcc_dF_restrict_equal;
phi_restrict = control.phi_restrict;
dphi_restrict = control.dphi_restrict;
ddphi_restrict = control.ddphi_restrict;
fcc_rose_equation = control.fcc_rose_equation;
fcc_100_tension = control.fcc_100_tension;
hcp_rose_equation = control.hcp_rose_equation;
bcc_rose_equation = control.bcc_rose_equation;
fcc_sf_relation = control.fcc_sf_relation;
fcc_surface_relation = control.fcc_surface_relation;

x0 = [F_data, phi_data];        % 优化初始解

% 关于FCC的一些基本参数
radius_fcc = alpha_fcc*a_fcc;               % FCC的半径
Aeq_phi_fcc = zeros(size(phi_data));        % FCC平衡点的phi矩阵
for i = 1:size(N_fcc, 2)
    Aeq_phi_fcc = Aeq_phi_fcc + 0.5*N_fcc(i)*cal_phi(radius_fcc(i), knots, order);
end
rho_fcc = sum(N_fcc .* interp1(rho(:,1), rho(:,2), radius_fcc));   % FCC平衡点的电荷密度
drho_fcc = interp1(drho_x, drho, radius_fcc);

% 给出F函数和phi函数的上限和下限
lb_F = -1*F_phi_range(1)*ones(size(F_data));
ub_F = F_phi_range(1)*ones(size(F_data));
lb_phi = -1*F_phi_range(2)*ones(size(phi_data));
ub_phi = F_phi_range(2)*ones(size(phi_data));
lb = [lb_F, lb_phi];     % lb的话，对于F是-1，对于phi是-100
ub = [ub_F, ub_phi];     % ub的话，对于F是1，对于phi是100

A = [];     % 初始化
b = [];

% 给出F函数的不等式约束,此处要求实测的C12-C44<0.5(C12_cons-C44_cons)的差距
% 很明显，这一条件仅适用于C12<C44的情况
if fcc_relation_C12_C44(1) == 1
    xi = rho_fcc;
    delta = [2, 2, 12, 8, 20];
    omega = a_fcc^3/4;
    V = sum(drho_fcc .* delta ./ alpha_fcc * a_fcc);
    A_temp(1, :) = [cal_ddF(xi), zeros(size(phi_data))];
    b_temp = fcc_relation_C12_C44(2) * (fcc_elastic(2)-fcc_elastic(3)) * omega/V^2;
    A = [A; A_temp];
    b = [b, b_temp];
end

% 保证F函数在[dF_restrict(3),dF_restrict(4)]的区间内是一阶导数小于0
if dF_restrict(1) == 1
    N_A = dF_restrict(2);
    xi = linspace(dF_restrict(3), dF_restrict(4), N_A);
    A_temp = zeros(N_A, size(x0, 2));
    for i = 1:N_A
        A_temp(i, :) = [cal_dF(xi(i)), zeros(size(phi_data))];
    end
    b_temp = 0 * ones(size(xi));        % 因为是小于0，,所以设计为0*ones
    A = [A; A_temp];
    b = [b, b_temp];
end

% 是否限制dF_rho_fcc值，限制的最大值，附近的采样点数
if fcc_dF_restrict_equal(1) == 1
    dF_cons = fcc_dF_restrict_equal(2);
    N_A = fcc_dF_restrict_equal(3);
    xi = linspace(rho_fcc-1, rho_fcc+1, N_A);
    A_temp = zeros(N_A, size(x0, 2));
    for i = 1:N_A
        A_temp(i, :) = [cal_dF(xi(i)), zeros(size(phi_data))];
    end
    b_temp = dF_cons * ones(size(xi));
    A = [A; A_temp];
    b = [b, b_temp];
end

% 给出phi函数在第一个knots处，最大的值
if phi_restrict(1) == 1
    A_temp = [zeros(size(F_data)), cal_phi(knots(1), knots, order)];
    b_temp = phi_restrict(2);
    F_phi_cons = F_phi_cons_record([0, 0], [1, knots(1)], [1, phi_restrict(2)], F_phi_cons);
    A = [A; A_temp];
    b = [b b_temp];
end


% 给出phi函数在不同位置处，一阶导数小于0
if dphi_restrict(1) == 1
    restrict_point = dphi_restrict(2:end);
    A_temp = zeros(size(restrict_point, 2), size(x0, 2));
    b_temp = zeros(1, size(restrict_point, 2));
    for i = 1:size(restrict_point, 2)
        A_temp(i, :) = [zeros(size(F_data)), cal_dphi(restrict_point(i), knots, order)];
        b_temp(i) = 0;
    end
    A = [A; A_temp];
    b = [b b_temp];
end

% 在一定的范围内，二阶导数不能超过一个限度，ddphi_restrict(2)
if ddphi_restrict(1) == 1
    N_A = ddphi_restrict(3);
    xi = linspace(ddphi_restrict(4), ddphi_restrict(5), N_A);
    A_temp = zeros(N_A, size(x0, 2));
    b_temp = ddphi_restrict(2) * ones(1, N_A);
    for i = 1:N_A
        A_temp(i, :) = [zeros(size(F_data)), cal_ddphi(xi(i), knots, order)];
    end
    A = [A; A_temp];
    b = [b b_temp];
    
    A_temp = zeros(N_A, size(x0, 2));
    b_temp = ddphi_restrict(2) * ones(1, N_A);
    for i = 1:N_A
        A_temp(i, :) = [zeros(size(F_data)), -1*cal_ddphi(xi(i), knots, order)];
    end
    A = [A; A_temp];
    b = [b b_temp];
end

% 给出等式约束条件，对FCC晶格常数的拟合，rose方程
if fcc_rose_equation(1) == 1
    N_equ = fcc_rose_equation(2);     % 选择几个点来做约束条件
    scan_span = [-fcc_rose_equation(3), fcc_rose_equation(3)];      % 扫描的区间
    r_wse = power(3/16/pi, 1/3) * a_fcc;         % 计算出等价体积
    K = (fcc_elastic(1) + 2*fcc_elastic(2))/3;   % 根据C11和C12计算出体积模量
    para_l = sqrt(sub_fcc/(12*pi*K*r_wse));     % 根据升华能(正值)和晶格常数，计算出中间参数l
    r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);
    a_norm = linspace(scan_span(1), scan_span(2), N_equ);
    energy_ref = sub_fcc * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);
    
    Aeq = zeros(N_equ, size(x0, 2));
    beq = zeros(1, N_equ);
    
    xi = linspace(-0.7, 0.7, N_equ);
    EPS_b = fcc_rose_equation(4)*(exp(6*abs(xi))/exp(6)*500+1);     % 给出不等式的容忍量
    for i = 1 : N_equ
        a_temp = power(16*pi/3, 1/3) * r(i);    % 适用于FCC晶格
        radius_this = alpha_fcc*a_temp;
        rho_this = sum(N_fcc.*interp1(rho(:,1), rho(:,2), radius_this));   % 计算出对应的电荷密度
        beq(i) = energy_ref(i);
        Aeq_phi = zeros(size(phi_data));
        for j = 1:size(radius_this, 2)
            Aeq_phi = Aeq_phi + 0.5*N_fcc(j)*cal_phi(radius_this(j), knots, order);
        end
        Aeq(i, :) = [cal_F(rho_this), Aeq_phi];

        F_coeff = [1, rho_this];
        phi_coeff = [0.5*N_fcc', radius_this'];
        b_coeff = [1, beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % 保存约束
        F_coeff = [-1, rho_this];
        phi_coeff = [-0.5*N_fcc', radius_this'];
        b_coeff = [1, -beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    end
    [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, fcc_rose_equation(5));   % 制作一个关于Rose曲线的不等式    
end

% 给出等式约束条件，对HCP晶格常数的拟合，rose方程。全程保持a变化，而c/a不变
if hcp_rose_equation(1) == 1
    N_equ = hcp_rose_equation(2);     % 选择几个点来做约束条件
    scan_span = [-hcp_rose_equation(3), hcp_rose_equation(3)];      % 扫描的区间
    r_wse = power(3*sqrt(3)/(32*pi)*(c_hcp/a_hcp), 1/3) * a_hcp;         % 计算出等价体积
    % K = (fcc_elastic(1) + 2*fcc_elastic(2))/3;   % 根据C11和C12计算出体积模量
    % 这里面需要使用FCC相的体积模量，因为以FCC为母相，为避免麻烦，直接采用FCC的体积模量
    
    para_l = sqrt(sub_hcp/(12*pi*K*r_wse));     % 根据升华能(正值)和晶格常数，计算出中间参数l
    r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);      % WS半径的扫描范围
    a_norm = linspace(scan_span(1), scan_span(2), N_equ);
    energy_ref = sub_hcp * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);    % 目标能量的参考值
    
    Aeq = zeros(N_equ, size(x0, 2));
    beq = zeros(1, N_equ);
    
    xi = linspace(-0.7, 0.7, N_equ);
    EPS_b = hcp_rose_equation(4)*(exp(6*abs(xi))/exp(6)*500+1);
    for i = 1 : N_equ
        a_temp = 1/power(3*sqrt(3)/(32*pi)*(c_hcp/a_hcp), 1/3) * r(i);    % 适用于HCP晶格
        c_temp = c_hcp/a_hcp * a_temp;
        [radius_this, N_this] = cal_hcp_radius(a_temp, c_temp, range(4));         % 获得新的HCP半径和原子数
        rho_this = sum(N_this.*interp1(rho(:,1), rho(:,2), radius_this));         % 计算出对应的电荷密度
        beq(i) = energy_ref(i);
        Aeq_phi = zeros(size(phi_data));
        for j = 1:size(radius_this, 2)
            Aeq_phi = Aeq_phi + 0.5*N_this(j)*cal_phi(radius_this(j), knots, order);
        end
        Aeq(i, :) = [cal_F(rho_this), Aeq_phi];
        
        F_coeff = [1, rho_this];
        phi_coeff = [0.5*N_this', radius_this'];
        b_coeff = [1, beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % 保存约束
        F_coeff = [-1, rho_this];
        phi_coeff = [-0.5*N_this', radius_this'];
        b_coeff = [1, -beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    end
    [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, hcp_rose_equation(5));   % 制作一个关于Rose曲线的不等式
end

% 给出等式约束条件，对BCC晶格常数的拟合，rose方程
if bcc_rose_equation(1) == 1
    N_equ = bcc_rose_equation(2);     % 选择几个点来做约束条件
    scan_span = [-bcc_rose_equation(3), bcc_rose_equation(3)];      % 扫描的区间
    r_wse = power(3/8/pi, 1/3) * a_bcc;         % 计算出等价体积
    % K = (fcc_elastic(1) + 2*fcc_elastic(2))/3;   % 根据C11和C12计算出体积模量
    % 这里面也是使用了FCC相的体积模量，这个程序仅适用于Fe50Mn30Cr10Co10的HEA合金
    
    para_l = sqrt(sub_bcc/(12*pi*K*r_wse));     % 根据升华能(正值)和晶格常数，计算出中间参数l
    r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);
    a_norm = linspace(scan_span(1), scan_span(2), N_equ);
    energy_ref = sub_bcc * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);
    
    Aeq = zeros(N_equ, size(x0, 2));
    beq = zeros(1, N_equ);
    
    xi = linspace(-0.7, 0.7, N_equ);
    EPS_b = bcc_rose_equation(4)*(exp(6*abs(xi))/exp(6)*500+1);
    for i = 1 : N_equ
        a_temp = power(8*pi/3, 1/3) * r(i);    % 适用于BCC晶格
        [radius_this, N_this] = cal_bcc_radius(a_temp, range(4));      % 计算出对应的几何因子
        rho_this = sum(N_this .* interp1(rho(:,1), rho(:,2), radius_this));   % 计算出对应的电荷密度
        beq(i) = energy_ref(i);
        Aeq_phi = zeros(size(phi_data));
        for j = 1:size(radius_this, 2)
            Aeq_phi = Aeq_phi + 0.5*N_this(j)*cal_phi(radius_this(j), knots, order);
        end
        Aeq(i, :) = [cal_F(rho_this), Aeq_phi];

        F_coeff = [1, rho_this];
        phi_coeff = [0.5*N_this', radius_this'];
        b_coeff = [1, beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % 保存约束
        F_coeff = [-1, rho_this];
        phi_coeff = [-0.5*N_this', radius_this'];
        b_coeff = [1, -beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    end
    [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, bcc_rose_equation(5));
end

% 给出针对于FCC晶格100方向拉伸的约束
if fcc_100_tension(1) == 1
    N_inequ = (size(fcc_100_tension, 2)-1)/2;
    index = 2;      % 数组的指针
    for i = 1 : N_inequ
        [A, b, F_phi_cons] = make_fcc_strain_100_inequal(A, b, rho, phi_data, Aeq_phi_fcc, rho_fcc, knots, order, a_fcc, fcc_100_tension(index:index+1), ...
            fcc_elastic(1), range(4), N_fcc, radius_fcc, F_phi_cons);
        index = index + 2;
    end
end

% 给出FCC晶体的不等式约束条件，要求广义平面缺陷能满足一般的曲线形式
% Egsf_0_1 < Egsf_0_2 < Egsf_0_3 < Egsf_0_4 < Egsf_0_5(不稳定层错能)
% Egsf_1_0(层错能) < Egsf_0_5(不稳定层错能)

if fcc_sf_relation(1) == 1
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [1, 2], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [2, 3], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [3, 4], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [4, 5], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [6, 5], F_phi_cons);
end

% 给出对于FCC表面能排序的约束，要求E111<E100<E110
if fcc_surface_relation(1) == 1
    [A, b, F_phi_cons] = make_fcc_surface_inequal(A, b, rho, phi_data, Aeq_phi_fcc, rho_fcc, knots, order, a_fcc, N_fcc_111, ...
        N_fcc_110, N_fcc_100, N_fcc, radius_fcc, F_phi_cons);
end

end
