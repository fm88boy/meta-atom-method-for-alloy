%% 读入fcc晶格的特殊构型与能量

function para = load_fcc_specific_config(para)
EPS = 1e-4;

if ~exist('specific_config.mat', 'file')
    error('缺少specific_config.mat文件！！！');
end
load specific_config            % 包含4个变量：target_lattice, cutoff, fcc_config, fcc_config_energy

if norm(para.cons_data.lattice - target_lattice) > EPS || abs(para.range(4)-cutoff) > EPS           % 如果特征晶格常数不匹配，则报错
    error('specific_config.mat的特征晶格常数不匹配！！！')
end

para.fcc_specific_config = fcc_config;
para.fcc_specific_energy = fcc_config_energy;
% config.geo(100, 1)
% config.F_geo(100, 1)
% energy(100, 1)
% config的数据结构：{[0.5, 6;
%                     1.5, 6;];     第一列几何参数，第二列原子个数
% energy的数据结构：[0.1; 0.2; 0.3; ......]       构型的相对能量

end