
function [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a, alpha_gsf, pair, F_phi_cons)
% 求出pair(1)<pair(2)的矩阵
alpha_1_1 = alpha_gsf{pair(1), 1};      % 层错面的第一层的几何因子
alpha_1_2 = alpha_gsf{pair(1), 2};      % 层错面的第二层的几何因子
alpha_2_1 = alpha_gsf{pair(2), 1};
alpha_2_2 = alpha_gsf{pair(2), 2};

rho_1_1 = sum(alpha_1_1(:, 2) .* interp1(rho(:,1), rho(:,2), alpha_1_1(:, 1)*a));   % 层错面第一层的电荷密度
rho_1_2 = sum(alpha_1_2(:, 2) .* interp1(rho(:,1), rho(:,2), alpha_1_2(:, 1)*a));
Aeq_phi_1_1 = zeros(size(phi_data));
for i = 1:size(alpha_1_1, 2)
    Aeq_phi_1_1 = Aeq_phi_1_1 + 0.5*alpha_1_1(i, 2) * cal_phi(alpha_1_1(i, 1)*a, knots, order);
end
Aeq_phi_1_2 = zeros(size(phi_data));
for i = 1:size(alpha_1_2, 2)
    Aeq_phi_1_2 = Aeq_phi_1_2 + 0.5*alpha_1_2(i, 2) * cal_phi(alpha_1_2(i, 1)*a, knots, order);
end
Aeq_phi_1 = Aeq_phi_1_1 + Aeq_phi_1_2;          % 因为有两层
Aeq_F_1 = cal_F(rho_1_1) + cal_F(rho_1_2);
Aeq_1 = [Aeq_F_1, Aeq_phi_1];

rho_2_1 = sum(alpha_2_1(:, 2) .* interp1(rho(:,1), rho(:,2), alpha_2_1(:, 1)*a));
rho_2_2 = sum(alpha_2_2(:, 2) .* interp1(rho(:,1), rho(:,2), alpha_2_2(:, 1)*a));
Aeq_phi_2_1 = zeros(size(phi_data));
for i = 1:size(alpha_2_1, 2)
    Aeq_phi_2_1 = Aeq_phi_2_1 + 0.5*alpha_2_1(i, 2)*cal_phi(alpha_2_1(i, 1)*a, knots, order);
end
Aeq_phi_2_2 = zeros(size(phi_data));
for i = 1:size(alpha_2_2, 2)
    Aeq_phi_2_2 = Aeq_phi_2_2 + 0.5*alpha_2_2(i, 2)*cal_phi(alpha_2_2(i, 1)*a, knots, order);
end
Aeq_phi_2 = Aeq_phi_2_1 + Aeq_phi_2_2;
Aeq_F_2 = cal_F(rho_2_1) + cal_F(rho_2_2);
Aeq_2 = [Aeq_F_2, Aeq_phi_2];

F_coeff = [1, rho_1_1;
           1, rho_1_2;
           -1, rho_2_1;
           -1, rho_2_2];
phi_coeff = [0.5*alpha_1_1(:, 2), alpha_1_1(:, 1)*a; 
             0.5*alpha_1_2(:, 2), alpha_1_2(:, 1)*a;
             -0.5*alpha_2_1(:, 2), alpha_2_1(:, 1)*a;
             -0.5*alpha_2_2(:, 2), alpha_2_2(:, 1)*a];
b_coeff = [1, 0];
F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % 保存约束

Aeq = Aeq_1 - Aeq_2;
beq = 0;
A = [A; Aeq];
b = [b beq];
end
