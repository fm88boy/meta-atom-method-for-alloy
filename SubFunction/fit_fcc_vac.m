%% 计算FCC的单原子空穴能

function [evac_fcc_cal] = fit_fcc_vac(F, rho, phi, N_fcc_vac, N_fcc, E_fcc, r_fcc)
rho_fcc_vac = zeros(1, 5);
phi_fcc_vac = zeros(1, 5);
rho_fcc_vac(1) = sum(N_fcc_vac(1, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_vac(1) = sum(N_fcc_vac(1, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
rho_fcc_vac(2) = sum(N_fcc_vac(2, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_vac(2) = sum(N_fcc_vac(2, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
rho_fcc_vac(3) = sum(N_fcc_vac(3, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_vac(3) = sum(N_fcc_vac(3, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
rho_fcc_vac(4) = sum(N_fcc_vac(4, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_vac(4) = sum(N_fcc_vac(4, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
rho_fcc_vac(5) = sum(N_fcc_vac(5, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_vac(5) = sum(N_fcc_vac(5, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
evac_fcc_cal = sum(N_fcc.*(interp1(F(:,1), F(:,2), rho_fcc_vac) + 0.5*phi_fcc_vac)) - sum(N_fcc)*E_fcc;
end