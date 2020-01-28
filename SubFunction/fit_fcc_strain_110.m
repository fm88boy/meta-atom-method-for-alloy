%% 计算fcc的110路径上的能量

function fcc_strain_110_cal = fit_fcc_strain_110(F, rho, phi, fcc_strain_110, a_fcc, Esub_fcc, cutoff)
% F，rho，phi：三个基础函数；fcc_strain_110：fcc晶格在110方向的应变和能量；
% a_fcc：FCC的晶格；Esub_fcc：FCC的升华能；cutoff：FCC势函数的截断半径

fcc_strain_110_cal = zeros(size(fcc_strain_110, 1), 2);
fcc_strain_110_cal(:, 1) = fcc_strain_110(:, 1);
alpha_fcc_strain_110{1} = [0.707106781187, 1.000000000000, 1.224744871392, 1.414213562373, 1.732050807569, 2.000000000000, 2.121320343560, 2.345207879912, 2.449489742783];
N_fcc_strain_110{1} = [2, 2, 4, 2, 4, 2, 6, 4, 4];
alpha_fcc_strain_110{2} = [0.612372435696, 1.172603939956, 1.541103500742, 1.837117307087, 2.318404623874];
N_fcc_strain_110{2} = [8, 8, 8, 16, 8];
alpha_fcc_strain_110{3} = [0.000000000000, 0.707106781187, 1.000000000000, 1.224744871392, 1.414213562373, 1.732050807569, 2.000000000000, 2.121320343560, 2.345207879912];
N_fcc_strain_110{3} = [2, 4, 4, 8, 4, 8, 4, 12, 8];
alpha_fcc_strain_110{4} = [0.612372435696, 1.172603939956, 1.541103500742, 1.837117307087];
N_fcc_strain_110{4} = [8, 8, 8, 16];
alpha_fcc_strain_110{5} = [0.000000000000, 0.707106781187, 1.000000000000, 1.224744871392, 1.414213562373, 1.732050807569, 2.000000000000];
N_fcc_strain_110{5} = [2, 4, 4, 8, 4, 8, 4];
alpha_fcc_strain_110{6} = [0.612372435696, 1.172603939956, 1.541103500742];
N_fcc_strain_110{6} = [8, 8, 8];
alpha_fcc_strain_110{7} = [0.000000000000, 0.707106781187, 1.000000000000, 1.224744871392];
N_fcc_strain_110{7} = [2, 4, 4, 8];
length = 41;

for i = 1:size(fcc_strain_110_cal, 1)
    [~, fcc_strain_110_cal(i, 2)] = loading_fcc_strain_110(F, rho, phi, fcc_strain_110_cal(i, 1), a_fcc, cutoff, ...
        alpha_fcc_strain_110, N_fcc_strain_110, length);
    fcc_strain_110_cal(i, 2) = fcc_strain_110_cal(i, 2) + Esub_fcc;
end
end
