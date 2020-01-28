%% ����FCC�ĵ��Գ���
function [C11_fcc_cal, C12_fcc_cal, C44_fcc_cal] = fit_fcc_elastic(ddF_fcc, dF_fcc, ddphi_fcc, dphi_fcc, ddrho_fcc, drho_fcc, ...
    alpha_fcc, a_fcc, critical, omega)
% ddF_fcc��dF_fcc��F����ƽ��λ�õĶ��׵�����һ�׵�����ddphi_fcc��dphi_fcc��phi����ƽ��λ�õĶ��׵�����һ�׵�����
% ddrho_fcc��drho_fcc��rho����ƽ��λ�õĶ��׵�����һ�׵�����
% alpha_fcc��FCC����ļ������ӣ�a_fcc��FCC�ľ���critical�����Գ����ļ������ӣ�omega������ԭ�����

B_fcc_1111 = 0.5*sum((ddphi_fcc-dphi_fcc./alpha_fcc/a_fcc) .* critical.alpha_fcc_1111 ./ alpha_fcc.^2 * a_fcc^2);
W_fcc_1111 = sum((ddrho_fcc - drho_fcc./(a_fcc*alpha_fcc)) .* critical.alpha_fcc_1111 ./ (alpha_fcc).^2 * a_fcc^2);
B_fcc_1122 = 0.5*sum((ddphi_fcc-dphi_fcc./alpha_fcc/a_fcc) .* critical.alpha_fcc_1122 ./ alpha_fcc.^2 * a_fcc^2);
W_fcc_1122 = sum((ddrho_fcc - drho_fcc./(a_fcc*alpha_fcc)) .* critical.alpha_fcc_1122 ./ (alpha_fcc).^2 * a_fcc^2);
V_fcc_11 = sum(drho_fcc .* critical.delta_fcc_11 ./ alpha_fcc * a_fcc);

C11_fcc_cal = (B_fcc_1111 + dF_fcc*W_fcc_1111 + ddF_fcc*V_fcc_11^2)/omega;
C12_fcc_cal = (B_fcc_1122 + dF_fcc*W_fcc_1122 + ddF_fcc*V_fcc_11^2)/omega;
C44_fcc_cal = (B_fcc_1122 + dF_fcc*W_fcc_1122)/omega;
end