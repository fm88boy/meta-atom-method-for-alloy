%% ����fcc��������⹹��������

function para = load_fcc_specific_config(para)
EPS = 1e-4;

if ~exist('specific_config.mat', 'file')
    error('ȱ��specific_config.mat�ļ�������');
end
load specific_config            % ����4��������target_lattice, cutoff, fcc_config, fcc_config_energy

if norm(para.cons_data.lattice - target_lattice) > EPS || abs(para.range(4)-cutoff) > EPS           % ���������������ƥ�䣬�򱨴�
    error('specific_config.mat��������������ƥ�䣡����')
end

para.fcc_specific_config = fcc_config;
para.fcc_specific_energy = fcc_config_energy;
% config.geo(100, 1)
% config.F_geo(100, 1)
% energy(100, 1)
% config�����ݽṹ��{[0.5, 6;
%                     1.5, 6;];     ��һ�м��β������ڶ���ԭ�Ӹ���
% energy�����ݽṹ��[0.1; 0.2; 0.3; ......]       ���͵��������

end