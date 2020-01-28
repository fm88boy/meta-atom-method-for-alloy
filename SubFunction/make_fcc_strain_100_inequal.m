
function [A, b, F_phi_cons] = make_fcc_strain_100_inequal(A, b, rho, phi_data, Aeq_phi_fcc, rho_fcc, knots, order, a_fcc, fcc_100_tension, ...
    C11, r_cut, N_fcc, radius_fcc, F_phi_cons)
% 计算几何因子
alpha_fcc_strain_100{1} = [0.707106781187, 1.000000000000, 1.414213562373, 1.581138830084, 2.000000000000, 2.121320343560, 2.236067977500];
N_fcc_strain_100{1} = [4, 4, 4, 8, 4, 4, 8];
alpha_fcc_strain_100{2} = [0.500000000000, 1.118033988750, 1.500000000000, 1.802775637732, 2.061552812809];
N_fcc_strain_100{2} = [8, 16, 8, 16, 16];
alpha_fcc_strain_100{3} = [0.000000000000, 0.707106781187, 1.000000000000, 1.414213562373, 1.581138830084, 2.000000000000, 2.121320343560, 2.236067977500];
N_fcc_strain_100{3} = [2, 8, 8, 8, 16, 8, 8, 16];
alpha_fcc_strain_100{4} = [0.500000000000, 1.118033988750, 1.500000000000, 1.802775637732];
N_fcc_strain_100{4} = [8, 16, 8, 16];
alpha_fcc_strain_100{5} = [0.000000000000, 0.707106781187, 1.000000000000, 1.414213562373];
N_fcc_strain_100{5} = [2, 8, 8, 8];
length = 28;

a_fcc_cal = a_fcc;                % 在横向上无应变
h_cal = 1/2*a_fcc * (1+fcc_100_tension(1));     % 这个变量需要根据方向而改变

radius = zeros(1, length);
N = zeros(1, length);
input_idx = 1;
for i = 1:size(alpha_fcc_strain_100, 2)
    this_length = size(alpha_fcc_strain_100{i}, 2);     % 该项长度
    radius(input_idx:input_idx+this_length-1) = sqrt((alpha_fcc_strain_100{i}*a_fcc_cal).^2 + ((i-1)*h_cal)^2);
    N(input_idx:input_idx+this_length-1) = N_fcc_strain_100{i};
    input_idx = input_idx + this_length;
end

idx = radius < r_cut;
radius = radius(idx);       % 对应的原子半径
N = N(idx);                 % 对应的原子个数

% 计算对应的能量
E_tot = 0.5*C11*fcc_100_tension(1)^2;
E_desire = fcc_100_tension(2)*E_tot*a_fcc^3/4;        % 乘以系数

% 给出约束条件
rho_this = sum(N .* interp1(rho(:,1), rho(:,2), radius));
Aeq_phi_this = zeros(size(phi_data));
for i = 1:size(radius, 2)
    Aeq_phi_this = Aeq_phi_this + 0.5 * N(i) * cal_phi(radius(i), knots, order);
end

F_coeff = [-1, rho_this;
           1, rho_fcc];
phi_coeff = [-0.5*N', radius';
             0.5*N_fcc', radius_fcc'];
b_coeff = [1, -1*E_desire];

F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % 保存约束

Aeq_phi = Aeq_phi_this - Aeq_phi_fcc;             % 相对fcc母相的能量
Aeq_F = cal_F(rho_this) - cal_F(rho_fcc);
Aeq = [Aeq_F, Aeq_phi];

Aeq = -1*Aeq;
beq = -1*E_desire;
A = [A; Aeq];
b = [b beq];
end
