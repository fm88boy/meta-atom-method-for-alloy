%% 计算phi函数的形状

function [diff_cal, max_phi_cal] = fit_phi_well(phi, phi_well)
% phi：phi函数；phi_well：控制函数形状的参数
% 
phi_cut_1 = phi(phi(:,1)>phi_well(1) & phi(:,1)<phi_well(2), :);    % 第一段区间
phi_cut_2 = phi(phi(:,1)>phi_well(2) & phi(:,1)<phi_well(3), :);    % 第二段区间

min_phi_1 = min(phi_cut_1(:, 2));
min_phi_2 = min(phi_cut_2(:, 2));
diff_cal = min_phi_2 - min_phi_1;

% 计算出震荡的幅度
phi_cut = phi(phi(:,1)>phi_well(5), :);      % 比第一近邻近0.1Ang
max_phi_cal = max(phi_cut(:, 2));   % 在最小值后面的最大值
end