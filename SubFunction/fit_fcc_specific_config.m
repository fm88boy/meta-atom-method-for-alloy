%% 计算特殊构型的能量

function fcc_specific_energy_cal = fit_fcc_specific_config(F, rho, phi, a_fcc, fcc_specific_config)
% F，rho，phi：三个基础函数；a_fcc：FCC的晶格；
% fcc_specific_config：fcc特殊构型的几何因子
num_config = size(fcc_specific_config.geo, 1);
fcc_specific_energy_cal = zeros(num_config, 1);

for i = 1:num_config
    fcc_specific_energy_cal(i) = cal_fcc_geo_energy(F, rho, phi, a_fcc, fcc_specific_config.geo{i}, fcc_specific_config.F_geo{i});
end
end