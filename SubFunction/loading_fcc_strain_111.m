%% 沿着111轴方向进行加载

function [e_expan, energy] = loading_fcc_strain_111(F, rho, phi, strain, a, r_cut, alpha_fcc_strain_111, N_fcc_strain_111, length)
% 计算了在一个晶格常数为a的FCC晶胞，在某个方向应变为strain时系统的能量，另外两个方向是自由放缩
bnd = [-1.5*abs(strain), 0.5*abs(strain)];    % 区间设计为-1.5*e_zz到0.5*e_zz，实际上，理论值应该是-0.5*e_zz
options = optimset('Display', 'off');
[e_expan, energy] = fminbnd(@(e_expan)expansion_fcc_111(F, rho, phi, a, e_expan, strain, r_cut, alpha_fcc_strain_111, N_fcc_strain_111, length), ...
    bnd(1), bnd(2), options);
end

function energy = expansion_fcc_111(F, rho, phi, a_fcc, e_expan, e_zz, r_cut, alpha_fcc_strain_111, N_fcc_strain_111, length)
a_fcc_cal = a_fcc * (1+e_expan);
h_cal = sqrt(3)/3*a_fcc * (1+e_zz);     % 这个变量需要根据方向而改变

radius = zeros(1, length);
N = zeros(1, length);
input_idx = 1;
for i = 1:size(alpha_fcc_strain_111, 2)
    this_length = size(alpha_fcc_strain_111{i}, 2);     % 该项长度
    radius(input_idx:input_idx+this_length-1) = sqrt((alpha_fcc_strain_111{i}*a_fcc_cal).^2 + ((i-1)*h_cal)^2);
    N(input_idx:input_idx+this_length-1) = N_fcc_strain_111{i};
    input_idx = input_idx + this_length;
end

idx = radius < r_cut;
radius = radius(idx);
N = N(idx);

phi_fcc = sum(N .* interp1(phi(:,1), phi(:,2), radius));
rho_fcc = sum(N .* interp1(rho(:,1), rho(:,2), radius));
energy = interp1(F(:, 1), F(:, 2), rho_fcc) + 0.5*phi_fcc;
end
