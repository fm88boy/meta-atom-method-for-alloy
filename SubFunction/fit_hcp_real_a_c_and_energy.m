%% 根据能量极小值的原理，找到hcp相的最低态的晶格常数与能量

function [a_hcp_cal, c_hcp_cal, hcp_energy_cal] = fit_hcp_real_a_c_and_energy(F, rho, phi, a_hcp, c_hcp, r_cut)
% F，rho，phi：三个基础函数；
% a_hcp, c_hcp：初始猜测的a_hcp, c_hcp常数；cutoff：势函数的截断半径
ratio = 0.08;
scan_step = 0.00001;                        % 扫描步长

lattice_init = [a_hcp, c_hcp];
lb = lattice_init * (1-ratio);
ub = lattice_init * (1+ratio);

options = optimset('Display', 'off', 'TolX', scan_step);
[lattice_hcp_cal, hcp_energy_cal] = fmincon(@(lattice_hcp)Inner_Lattice_Consistent_HCP(F, rho, phi, lattice_hcp, r_cut), ...
    lattice_init, [], [], [], [], lb, ub, [], options);
a_hcp_cal = lattice_hcp_cal(1);
c_hcp_cal = lattice_hcp_cal(2);
end

function energy = Inner_Lattice_Consistent_HCP(F, rho, phi, lattice_hcp, r_cut)
[r_hcp, N_hcp] = cal_hcp_radius(lattice_hcp(1), lattice_hcp(2), r_cut);   % 计算出HCP相的原子半径和原子个数
phi_hcp = sum(N_hcp .* interp1(phi(:,1), phi(:,2), r_hcp));
rho_hcp = sum(N_hcp .* interp1(rho(:,1), rho(:,2), r_hcp));
energy = interp1(F(:, 1), F(:, 2), rho_hcp) + 0.5*phi_hcp;
end
