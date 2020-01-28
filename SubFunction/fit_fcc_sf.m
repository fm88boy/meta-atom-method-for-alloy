%% 计算FCC的层错能与不稳定层错能

function [esf_fcc_cal, eusf_fcc_cal] = fit_fcc_sf(F, rho, phi, alpha_fcc_sf, N_fcc_sf, alpha_fcc_usf1, N_fcc_usf1, ...
    alpha_fcc_usf2, N_fcc_usf2, a_fcc, E_fcc)

rho_fcc_sf = sum(N_fcc_sf .* interp1(rho(:,1), rho(:,2), alpha_fcc_sf*a_fcc));
phi_fcc_sf = sum(N_fcc_sf .* interp1(phi(:,1), phi(:,2), alpha_fcc_sf*a_fcc));
E_fcc_sf = interp1(F(:,1), F(:,2), rho_fcc_sf) + 0.5*phi_fcc_sf;
esf_fcc_cal = 16*sqrt(3)/3*(E_fcc_sf - E_fcc)/a_fcc^2;

rho_fcc_usf1 = sum(N_fcc_usf1 .* interp1(rho(:,1), rho(:,2), alpha_fcc_usf1 * a_fcc));
phi_fcc_usf1 = sum(N_fcc_usf1 .* interp1(phi(:,1), phi(:,2), alpha_fcc_usf1 * a_fcc));
E_fcc_usf1 = interp1(F(:, 1), F(:, 2), rho_fcc_usf1) + 0.5 * phi_fcc_usf1;
rho_fcc_usf2 = sum(N_fcc_usf2 .* interp1(rho(:,1), rho(:,2), alpha_fcc_usf2 * a_fcc));
phi_fcc_usf2 = sum(N_fcc_usf2 .* interp1(phi(:,1), phi(:,2), alpha_fcc_usf2 * a_fcc));
E_fcc_usf2 = interp1(F(:, 1), F(:, 2), rho_fcc_usf2) + 0.5 * phi_fcc_usf2;
eusf_fcc_cal = 8*sqrt(3)/3*(E_fcc_usf1 + E_fcc_usf2 - 2*E_fcc)/a_fcc^2;
end