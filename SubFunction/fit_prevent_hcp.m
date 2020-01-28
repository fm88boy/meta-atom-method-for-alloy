%% 计算阻碍HCP相变的能量

function hcp_energy = fit_prevent_hcp(F, rho, phi, z_range, a_fcc, Esub_fcc, cutoff)
a_hcp = a_fcc/sqrt(2);
c_hcp = 2*sqrt(3)/3*a_fcc;

[~, hcp_energy] = loading_hcp(F, rho, phi, z_range, a_hcp, c_hcp, cutoff);
hcp_energy = hcp_energy + Esub_fcc;
end

