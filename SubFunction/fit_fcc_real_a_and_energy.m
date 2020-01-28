%% 根据能量极小值的原理，找到fcc相的最低态的晶格常数与能量

function [a_fcc_cal, fcc_energy_cal] = fit_fcc_real_a_and_energy(F, rho, phi, a_fcc, r_cut)
% F，rho，phi：三个基础函数；
% a_fcc：初始猜测的a_fcc常数；cutoff：势函数的截断半径
alpha_fcc = [0.707106781187, 1.000000000000, 1.224744871392, 1.414213562373, 1.581138830084, 1.732050807569, 1.870828693387, ...
    2.000000000000, 2.121320343560, 2.236067977500, 2.345207879912, 2.449489742783];
N_fcc = [12, 6, 24, 12, 24, 8, 48, 6, 36, 24, 24, 24];

ratio = 0.08;
a_range = a_fcc * [1-ratio, 1+ratio];       % 扫描范围
scan_step = 0.00001;                        % 扫描步长

options = optimset('Display', 'off', 'TolX', scan_step);
[a_fcc_cal, fcc_energy_cal] = fminbnd(@(a_fcc)Inner_Lattice_Consistent_FCC(F, rho, phi, alpha_fcc, N_fcc, a_fcc, r_cut), a_range(1), a_range(2), options);
end

function energy = Inner_Lattice_Consistent_FCC(F, rho, phi, alpha_fcc, N_fcc, a_fcc, r_cut)
radius = alpha_fcc * a_fcc;
N = N_fcc;

idx = radius < r_cut;
radius = radius(idx);
N = N(idx);

phi_fcc = sum(N .* interp1(phi(:,1), phi(:,2), radius));
rho_fcc = sum(N .* interp1(rho(:,1), rho(:,2), radius));
energy = interp1(F(:, 1), F(:, 2), rho_fcc) + 0.5*phi_fcc;
end
