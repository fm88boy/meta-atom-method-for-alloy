%% 计算FCC的表面能，分别为111,110和100面

function [e111_fcc_cal, e110_fcc_cal, e100_fcc_cal] = fit_fcc_surface(F, rho, phi, N_fcc_111, N_fcc_110, N_fcc_100, r_fcc, a_fcc, E_fcc)
rho_fcc_111_1 = sum(N_fcc_111(1, :) .* interp1(rho(:, 1), rho(:, 2), r_fcc));
phi_fcc_111_1 = sum(N_fcc_111(1, :) .* interp1(phi(:, 1), phi(:, 2), r_fcc));
E_fcc_111_1 = interp1(F(:, 1), F(:, 2), rho_fcc_111_1) + 0.5 * phi_fcc_111_1;
rho_fcc_111_2 = sum(N_fcc_111(2, :) .* interp1(rho(:, 1), rho(:, 2), r_fcc));
phi_fcc_111_2 = sum(N_fcc_111(2, :) .* interp1(phi(:, 1), phi(:, 2), r_fcc));
E_fcc_111_2 = interp1(F(:, 1), F(:, 2), rho_fcc_111_2) + 0.5 * phi_fcc_111_2;
e111_fcc_cal = 4*sqrt(3)/3*(E_fcc_111_1 + E_fcc_111_2 - 2*E_fcc)/a_fcc^2;

rho_fcc_110_1 = sum(N_fcc_110(1, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_110_1 = sum(N_fcc_110(1, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_110_1 = interp1(F(:, 1), F(:, 2), rho_fcc_110_1) + 0.5 * phi_fcc_110_1;
rho_fcc_110_2 = sum(N_fcc_110(2, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_110_2 = sum(N_fcc_110(2, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_110_2 = interp1(F(:, 1), F(:, 2), rho_fcc_110_2) + 0.5 * phi_fcc_110_2;
rho_fcc_110_3 = sum(N_fcc_110(3, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_110_3 = sum(N_fcc_110(3, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_110_3 = interp1(F(:, 1), F(:, 2), rho_fcc_110_3) + 0.5 * phi_fcc_110_3;
rho_fcc_110_4 = sum(N_fcc_110(4, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_110_4 = sum(N_fcc_110(4, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_110_4 = interp1(F(:, 1), F(:, 2), rho_fcc_110_4) + 0.5 * phi_fcc_110_4;
e110_fcc_cal = sqrt(2)*(E_fcc_110_1 + E_fcc_110_2 + E_fcc_110_3 + E_fcc_110_4 - 4*E_fcc)/a_fcc^2;

rho_fcc_100_1 = sum(N_fcc_100(1, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_100_1 = sum(N_fcc_100(1, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_100_1 = interp1(F(:, 1), F(:, 2), rho_fcc_100_1) + 0.5 * phi_fcc_100_1;
rho_fcc_100_2 = sum(N_fcc_100(2, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_100_2 = sum(N_fcc_100(2, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_100_2 = interp1(F(:, 1), F(:, 2), rho_fcc_100_2) + 0.5 * phi_fcc_100_2;
rho_fcc_100_3 = sum(N_fcc_100(3, :) .* interp1(rho(:,1), rho(:,2), r_fcc));
phi_fcc_100_3 = sum(N_fcc_100(3, :) .* interp1(phi(:,1), phi(:,2), r_fcc));
E_fcc_100_3 = interp1(F(:, 1), F(:, 2), rho_fcc_100_3) + 0.5 * phi_fcc_100_3;
e100_fcc_cal = 2*(E_fcc_100_1 + E_fcc_100_2 + E_fcc_100_3 - 3*E_fcc)/a_fcc^2;
end