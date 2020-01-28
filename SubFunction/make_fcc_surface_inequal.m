%% 制作关于FCC表面能的约束条件

function [A, b, F_phi_cons] = make_fcc_surface_inequal(A, b, rho, phi_data, Aeq_phi_fcc, rho_fcc, knots, order, a_fcc, N_fcc_111, ...
    N_fcc_110, N_fcc_100, N_fcc, radius_fcc, F_phi_cons)

[Aeq_111, F_coeff_111, phi_coeff_111] = cal_surface(4*sqrt(3)/3/a_fcc^2, phi_data, radius_fcc, Aeq_phi_fcc, rho_fcc, N_fcc_111, rho, knots, order, N_fcc);
[Aeq_110, F_coeff_110, phi_coeff_110] = cal_surface(sqrt(2)/a_fcc^2, phi_data, radius_fcc, Aeq_phi_fcc, rho_fcc, N_fcc_110, rho, knots, order, N_fcc);
[Aeq_100, F_coeff_100, phi_coeff_100] = cal_surface(2/a_fcc^2, phi_data, radius_fcc, Aeq_phi_fcc, rho_fcc, N_fcc_100, rho, knots, order, N_fcc);

F_coeff_sub = [F_coeff_111; -1*F_coeff_100(:, 1), F_coeff_100(:, 2)];
phi_coeff_sub = [phi_coeff_111; -1*phi_coeff_100(:, 1), phi_coeff_100(:, 2)];
b_coeff_sub = [1, 0];
F_phi_cons = F_phi_cons_record(F_coeff_sub, phi_coeff_sub, b_coeff_sub, F_phi_cons);    % 保存约束

A_temp = Aeq_111 - Aeq_100;
b_temp = 0;
A = [A; A_temp];
b = [b b_temp];

F_coeff_sub = [F_coeff_100; -1*F_coeff_110(:, 1), F_coeff_110(:, 2)];
phi_coeff_sub = [phi_coeff_100; -1*phi_coeff_110(:, 1), phi_coeff_110(:, 2)];
b_coeff_sub = [1, 0];
F_phi_cons = F_phi_cons_record(F_coeff_sub, phi_coeff_sub, b_coeff_sub, F_phi_cons);    % 保存约束

A_temp = Aeq_100 - Aeq_110;
b_temp = 0;
A = [A; A_temp];
b = [b b_temp];
end

%% 根据表面的几何因子，计算表面能的矩阵
function [Aeq_surface, F_coeff, phi_coeff] = cal_surface(coeff, phi_data, radius, Aequ_phi_lattice, rho_lattice, N_surface, rho, knots, order, N_lattice)
% coeff：量纲转化系数，phi_data：phi函数的参数，a_lattice：晶格常数，alpha_lattice：完美晶格的几何因子
% Aequ_phi_lattice：计算平衡态能量的矩阵，rho_lattice：完美晶格的电子密度，N_surface：表面的几何因子
% rho：电子密度，knots：phi函数的节点，order：phi函数的阶次
F_coeff = zeros(100, 2);
phi_coeff = zeros(500, 2);
F_idx = 1;
phi_idx = 1;

Aeq_phi = zeros(size(phi_data));
Aeq_F = 0;
for i = 1:size(N_surface, 1)    % 根据表面几何因子的层数来区分
    N_layer = N_surface(i, :);
    rho_layer = sum(N_layer .* interp1(rho(:,1), rho(:,2), radius));
    Aeq_phi_layer = zeros(size(phi_data));
    for j = 1:size(N_layer, 2)
        Aeq_phi_layer = Aeq_phi_layer + 0.5*N_layer(j)*cal_phi(radius(j), knots, order);
    end
    Aeq_phi = Aeq_phi + Aeq_phi_layer - Aequ_phi_lattice;
    Aeq_F = Aeq_F + cal_F(rho_layer) - cal_F(rho_lattice);
    F_coeff(F_idx:F_idx+1, :) = [1, rho_layer; -1, rho_lattice];
    phi_coeff(phi_idx:phi_idx+2*size(radius, 2)-1, :) = [0.5*N_layer', radius'; -0.5*N_lattice', radius'];
    F_idx = F_idx + 2;
    phi_idx = phi_idx + 2*size(radius, 2);
end
F_coeff(F_idx:end, :) = [];         % 将无用的数据删除
phi_coeff(phi_idx:end, :) = [];

Aeq_surface = coeff * [Aeq_F, Aeq_phi];     % 在系数上，乘以面积的折算因子
F_coeff(:, 1) = coeff * F_coeff(:, 1);
phi_coeff(:, 1) = coeff * phi_coeff(:, 1);
end