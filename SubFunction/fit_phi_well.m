%% ����phi��������״

function [diff_cal, max_phi_cal] = fit_phi_well(phi, phi_well)
% phi��phi������phi_well�����ƺ�����״�Ĳ���
% 
phi_cut_1 = phi(phi(:,1)>phi_well(1) & phi(:,1)<phi_well(2), :);    % ��һ������
phi_cut_2 = phi(phi(:,1)>phi_well(2) & phi(:,1)<phi_well(3), :);    % �ڶ�������

min_phi_1 = min(phi_cut_1(:, 2));
min_phi_2 = min(phi_cut_2(:, 2));
diff_cal = min_phi_2 - min_phi_1;

% ������𵴵ķ���
phi_cut = phi(phi(:,1)>phi_well(5), :);      % �ȵ�һ���ڽ�0.1Ang
max_phi_cal = max(phi_cut(:, 2));   % ����Сֵ��������ֵ
end