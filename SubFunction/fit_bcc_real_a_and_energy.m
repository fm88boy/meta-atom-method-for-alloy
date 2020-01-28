%% 根据能量极小值的原理，找到bcc相的最低态的晶格常数与能量

function [a_bcc_cal, bcc_energy_cal] = fit_bcc_real_a_and_energy(F, rho, phi, a_bcc, r_cut)
% F，rho，phi：三个基础函数；
% a_bcc：初始猜测的a_bcc常数；cutoff：势函数的截断半径
alpha_bcc = [0.866025403784, 1.000000000000, 1.414213562373, 1.658312395178, 1.732050807569, 2.000000000000, ...
    2.179449471770, 2.236067977500, 2.449489742783, 2.598076211353, 2.828427124746, 2.958039891550, 3.000000000000, 3.162277660168];
N_bcc = [8, 6, 12, 24, 8, 6, 24, 24, 24, 32, 12, 48, 30, 24];

ratio = 0.08;
a_range = a_bcc * [1-ratio, 1+ratio];       % 扫描范围
scan_step = 0.00001;                        % 扫描步长

options = optimset('Display', 'off', 'TolX', scan_step);
[a_bcc_cal, bcc_energy_cal] = fminbnd(@(a_bcc)Inner_Lattice_Consistent_BCC(F, rho, phi, alpha_bcc, N_bcc, a_bcc, r_cut), a_range(1), a_range(2), options);
end

function energy = Inner_Lattice_Consistent_BCC(F, rho, phi, alpha_bcc, N_bcc, a_bcc, r_cut)
radius = alpha_bcc * a_bcc;
N = N_bcc;

idx = radius < r_cut;
radius = radius(idx);
N = N(idx);

phi_bcc = sum(N .* interp1(phi(:,1), phi(:,2), radius));
rho_bcc = sum(N .* interp1(rho(:,1), rho(:,2), radius));
energy = interp1(F(:, 1), F(:, 2), rho_bcc) + 0.5*phi_bcc;
end
