%% 计算FCC相在拉伸时的能量

function stretch_energy = fit_maintain_fcc_energy(F, rho, phi, maintain_fcc_energy, a_fcc, Esub_fcc, cutoff)
% F,rho,phi：三个基础函数；maintain_fcc_energy：第一列应变，第二列能量
% a_fcc，Esub_fcc：fcc相的晶格常数，升华能；cutoff：势函数的截断半径
stretch_energy = zeros(size(maintain_fcc_energy, 1), 2);
for i = 1:size(maintain_fcc_energy, 1)
    [~, fcc_energy] = loading_fcc(F, rho, phi, maintain_fcc_energy(i, 1), a_fcc, cutoff);
    stretch_energy(i, :) = [maintain_fcc_energy(i, 1), fcc_energy+Esub_fcc];
end
end

