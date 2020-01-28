%% ���hcp����ʵrose����

function [a_real, energy_ref, energy] = fit_hcp_rose(F, rho, phi, N_equ, scan_span, a_hcp, c_hcp, Esub_hcp, K, r_cut)
r_wse = power(3*sqrt(3)/(32*pi)*(c_hcp/a_hcp), 1/3) * a_hcp;         % ������ȼ����
para_l = sqrt(Esub_hcp/(12*pi*K*r_wse));     % ����������(��ֵ)�;�������������м����l
r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);      % WS�뾶��ɨ�跶Χ

a_norm = linspace(scan_span(1), scan_span(2), N_equ);
energy_ref = Esub_hcp * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);    % Ŀ�������Ĳο�ֵ
a_real = 1/power(3*sqrt(3)/(32*pi)*(c_hcp/a_hcp), 1/3) * r;

energy = zeros(size(energy_ref));
for i = 1:size(energy, 2)
    a_temp = a_real(i);            % ������HCP����
    c_temp = c_hcp/a_hcp * a_temp;
    [radius_temp, N_temp] = cal_hcp_radius(a_temp, c_temp, r_cut);         % ����µ�HCP�뾶��ԭ����
    rho_temp = sum(N_temp .* interp1(rho(:,1), rho(:,2), radius_temp));         % �������Ӧ�ĵ���ܶ�
    energy(i) = interp1(F(:,1), F(:,2), rho_temp) + 0.5 * sum(N_temp.*interp1(phi(:,1), phi(:,2), radius_temp));
end
end
