%% 计算fcc的100路径（Bain路径）上的能量

function fcc_strain_100_cal = fit_fcc_strain_100(F, rho, phi, fcc_strain_100, a_fcc, Esub_fcc, cutoff)
% F，rho，phi：三个基础函数；fcc_strain_100：fcc晶格在100方向的应变和能量；
% a_fcc：FCC的晶格；Esub_fcc：FCC的升华能；cutoff：FCC势函数的截断半径

fcc_strain_100_cal = zeros(size(fcc_strain_100, 1), 2);
fcc_strain_100_cal(:, 1) = fcc_strain_100(:, 1);
alpha_fcc_strain_100{1} = [0.707106781187, 1.000000000000, 1.414213562373, 1.581138830084, 2.000000000000, 2.121320343560, 2.236067977500];
N_fcc_strain_100{1} = [4, 4, 4, 8, 4, 4, 8];
alpha_fcc_strain_100{2} = [0.500000000000, 1.118033988750, 1.500000000000, 1.802775637732, 2.061552812809];
N_fcc_strain_100{2} = [8, 16, 8, 16, 16];
alpha_fcc_strain_100{3} = [0.000000000000, 0.707106781187, 1.000000000000, 1.414213562373, 1.581138830084, 2.000000000000, 2.121320343560, 2.236067977500];
N_fcc_strain_100{3} = [2, 8, 8, 8, 16, 8, 8, 16];
alpha_fcc_strain_100{4} = [0.500000000000, 1.118033988750, 1.500000000000, 1.802775637732];
N_fcc_strain_100{4} = [8, 16, 8, 16];
alpha_fcc_strain_100{5} = [0.000000000000, 0.707106781187, 1.000000000000, 1.414213562373];
N_fcc_strain_100{5} = [2, 8, 8, 8];
length = 28;

for i = 1:size(fcc_strain_100_cal, 1)
    [~, fcc_strain_100_cal(i, 2)] = loading_fcc_strain_100(F, rho, phi, fcc_strain_100_cal(i, 1), a_fcc, cutoff, ...
        alpha_fcc_strain_100, N_fcc_strain_100, length);
    fcc_strain_100_cal(i, 2) = fcc_strain_100_cal(i, 2) + Esub_fcc;
end
end
