%% 计算fcc的111路径上的能量

function fcc_strain_111_cal = fit_fcc_strain_111(F, rho, phi, fcc_strain_111, a_fcc, Esub_fcc, cutoff)
% F，rho，phi：三个基础函数；fcc_strain_111：fcc晶格在111方向的应变和能量；
% a_fcc：FCC的晶格；Esub_fcc：FCC的升华能；cutoff：FCC势函数的截断半径

fcc_strain_111_cal = zeros(size(fcc_strain_111, 1), 2);
fcc_strain_111_cal(:, 1) = fcc_strain_111(:, 1);
alpha_fcc_strain_111{1} = [0.707106781187, 1.224744871392, 1.414213562373, 1.870828693387, 2.121320343560, 2.449489742783];
N_fcc_strain_111{1} = [6, 6, 6, 12, 6, 6];
alpha_fcc_strain_111{2} = [0.408248290464, 0.816496580928, 1.080123449735, 1.471960144388, 1.632993161855, 1.779513042005, 2.041241452319, 2.160246899469, 2.273030282831];
N_fcc_strain_111{2} = [6, 6, 12, 12, 6, 12, 6, 12, 12];
alpha_fcc_strain_111{3} = [0.408248290464, 0.816496580928, 1.080123449735, 1.471960144388, 1.632993161855, 1.779513042005, 2.041241452319, 2.160246899469];
N_fcc_strain_111{3} = [6, 6, 12, 12, 6, 12, 6, 12];
alpha_fcc_strain_111{4} = [0.000000000000, 0.707106781187, 1.224744871392, 1.414213562373];
N_fcc_strain_111{4} = [2, 12, 12, 12];
alpha_fcc_strain_111{5} = [0.408248290464, 0.816496580928];
N_fcc_strain_111{5} = [6, 6];
length = 29;

for i = 1:size(fcc_strain_111_cal, 1)
    [~, fcc_strain_111_cal(i, 2)] = loading_fcc_strain_111(F, rho, phi, fcc_strain_111_cal(i, 1), a_fcc, cutoff, ...
        alpha_fcc_strain_111, N_fcc_strain_111, length);
    fcc_strain_111_cal(i, 2) = fcc_strain_111_cal(i, 2) + Esub_fcc;
end
end
