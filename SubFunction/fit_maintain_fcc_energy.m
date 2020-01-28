%% ����FCC��������ʱ������

function stretch_energy = fit_maintain_fcc_energy(F, rho, phi, maintain_fcc_energy, a_fcc, Esub_fcc, cutoff)
% F,rho,phi����������������maintain_fcc_energy����һ��Ӧ�䣬�ڶ�������
% a_fcc��Esub_fcc��fcc��ľ������������ܣ�cutoff���ƺ����Ľضϰ뾶
stretch_energy = zeros(size(maintain_fcc_energy, 1), 2);
for i = 1:size(maintain_fcc_energy, 1)
    [~, fcc_energy] = loading_fcc(F, rho, phi, maintain_fcc_energy(i, 1), a_fcc, cutoff);
    stretch_energy(i, :) = [maintain_fcc_energy(i, 1), fcc_energy+Esub_fcc];
end
end

