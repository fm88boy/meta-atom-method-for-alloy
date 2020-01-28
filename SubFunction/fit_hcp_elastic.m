function [C11_hcp_cal, C12_hcp_cal, C13_hcp_cal, C33_hcp_cal, C44_hcp_cal] = fit_hcp_elastic(ddF_hcp, dF_hcp, ddphi_hcp, dphi_hcp, ...
    ddrho_hcp, drho_hcp, alpha_hcp, a_fcc, critical, omega)
error('该函数依赖于FCC相，并不是真实的HCP弹性常数程序')
B_hcp_1111 = 0.5*sum((ddphi_hcp-dphi_hcp./alpha_hcp/a_fcc) .* critical.alpha_hcp_1111 ./ alpha_hcp.^2 * a_fcc^2);
W_hcp_1111 = sum((ddrho_hcp - drho_hcp./(a_fcc*alpha_hcp)) .* critical.alpha_hcp_1111 ./ (alpha_hcp).^2 * a_fcc^2);
B_hcp_1122 = 0.5*sum((ddphi_hcp-dphi_hcp./alpha_hcp/a_fcc) .* critical.alpha_hcp_1122 ./ alpha_hcp.^2 * a_fcc^2);
W_hcp_1122 = sum((ddrho_hcp - drho_hcp./(a_fcc*alpha_hcp)) .* critical.alpha_hcp_1122 ./ (alpha_hcp).^2 * a_fcc^2);
B_hcp_1133 = 0.5*sum((ddphi_hcp-dphi_hcp./alpha_hcp/a_fcc) .* critical.alpha_hcp_1133 ./ alpha_hcp.^2 * a_fcc^2);
W_hcp_1133 = sum((ddrho_hcp - drho_hcp./(a_fcc*alpha_hcp)) .* critical.alpha_hcp_1133 ./ (alpha_hcp).^2 * a_fcc^2);
B_hcp_3333 = 0.5*sum((ddphi_hcp-dphi_hcp./alpha_hcp/a_fcc) .* critical.alpha_hcp_3333 ./ alpha_hcp.^2 * a_fcc^2);
W_hcp_3333 = sum((ddrho_hcp - drho_hcp./(a_fcc*alpha_hcp)) .* critical.alpha_hcp_3333 ./ (alpha_hcp).^2 * a_fcc^2);
V_hcp_11 = sum(drho_hcp .* critical.delta_hcp_11 ./ alpha_hcp * a_fcc);
V_hcp_33 = sum(drho_hcp .* critical.delta_hcp_33 ./ alpha_hcp * a_fcc);
C11_hcp_cal = (B_hcp_1111 + dF_hcp*W_hcp_1111 + ddF_hcp*V_hcp_11^2)/omega;
C12_hcp_cal = (B_hcp_1122 + dF_hcp*W_hcp_1122 + ddF_hcp*V_hcp_11^2)/omega;
C13_hcp_cal = (B_hcp_1133 + dF_hcp*W_hcp_1133 + ddF_hcp*V_hcp_11*V_hcp_33)/omega;
C33_hcp_cal = (B_hcp_3333 + dF_hcp*W_hcp_3333 + ddF_hcp*V_hcp_33^2)/omega;
C44_hcp_cal = (B_hcp_1133 + dF_hcp*W_hcp_1133)/omega;
end