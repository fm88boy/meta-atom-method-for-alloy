%% 计算BCC的Rose曲线

function [a_real, energy_ref, energy] = fit_bcc_rose(F, rho, phi, N_equ, scan_span, a_bcc, esub_bcc, K_bcc, r_cut)
r_wse = power(3/8/pi, 1/3) * a_bcc;         % 计算出等价体积
l = sqrt(esub_bcc/(12*pi*K_bcc*r_wse));
r = linspace(r_wse+scan_span(1)*l, r_wse+scan_span(2)*l, N_equ);

a_norm = linspace(scan_span(1), scan_span(2), N_equ);
energy_ref = esub_bcc * (-1 - a_norm - 0.05*a_norm.^3) .* exp(-a_norm);
a_real = power(8*pi/3, 1/3) * (a_norm*l + r_wse);

energy = zeros(size(r));
for i = 1:size(r, 2)
    a_temp = power(8*pi/3, 1/3) * r(i);    % 适用于BCC晶格
    [radius_this, N_this] = cal_bcc_radius(a_temp, r_cut);      % 计算出对应的几何因子
    rho_temp = sum(N_this .* interp1(rho(:,1), rho(:,2), radius_this));   % 计算出对应的电荷密度
    energy(i) = interp1(F(:,1), F(:,2), rho_temp) + 0.5*sum(N_this.*interp1(phi(:,1), phi(:,2), radius_this));
end
end