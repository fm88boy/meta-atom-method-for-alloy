%% ¼ÆËãFCCµÄRoseÇúÏß

function [a_real, energy_ref, energy] = fit_fcc_rose(F, rho, phi, N_equ, scan_span, a_fcc, esub_fcc, C11_fcc, C12_fcc, alpha_fcc, N_fcc)
K_fcc = (C11_fcc + 2*C12_fcc)/3;
r_wse = power(3/16/pi, 1/3) * a_fcc;
l = sqrt(esub_fcc/(12*pi * K_fcc * r_wse));
r = linspace(r_wse+scan_span(1)*l, r_wse+scan_span(2)*l, N_equ);

a_norm = linspace(scan_span(1), scan_span(2), N_equ);
energy_ref = esub_fcc * (-1 - a_norm - 0.05*a_norm.^3) .* exp(-a_norm);
a_real = power(16*pi/3, 1/3) * (a_norm*l + r_wse);

energy = zeros(size(r));
for i = 1:size(r, 2)
    a_temp = power(16*pi/3, 1/3) * r(i);
    rho_temp = sum(N_fcc.*interp1(rho(:,1), rho(:,2), alpha_fcc*a_temp));
    energy(i) = interp1(F(:,1), F(:,2), rho_temp) + 0.5*sum(N_fcc.*interp1(phi(:,1), phi(:,2), alpha_fcc*a_temp));
end
end