%% 计算FCC相的声子频率

function frequency_fcc_cal = fit_fcc_frequency(D, ddphi_fcc, dphi_fcc, ddF_fcc, dF_fcc, ddrho_fcc, drho_fcc, alpha_fcc, r_fcc)
% D：需要优化声子点的动态矩阵
% ddF_fcc，dF_fcc：F函数平衡位置的二阶导数和一阶导数；ddphi_fcc，dphi_fcc：phi函数平衡位置的二阶导数和一阶导数；
% alpha_fcc：FCC晶格的几何因子；r_fcc：FCC晶格的平衡位置

Dyn = zeros(size(D, 1), 6);
P = zeros(size(D, 1), 3);

for k = 1:size(D, 1)
    for i = 1:size(alpha_fcc, 2)
        temp = D{k}(i, 1:6);
        Dyn(k, 1:6) = Dyn(k, 1:6) + ddphi_fcc(i) * temp - 1/r_fcc(i) * dphi_fcc(i) * temp ...
            + 2*dF_fcc*(ddrho_fcc(i) * temp - 1/r_fcc(i) * drho_fcc(i) * temp);
        Dyn(k, 1:3) = Dyn(k, 1:3) + 1/r_fcc(i) * dphi_fcc(i) * D{k}(i, 7) + 2*dF_fcc*(1/r_fcc(i) * drho_fcc(i) * D{k}(i, 7));
        P(k,:) = P(k,:) + drho_fcc(i) * D{k}(i, 8:10);
    end
    temp = P(k, :)' * P(k, :);
    Dyn(k, :) = Dyn(k, :) - ddF_fcc * [temp(1,1), temp(2,2), temp(3,3), temp(1,2), temp(2,3), temp(1,3)];
end

Dyn(:,1) = Dyn(:,1) - Dyn(1,1);
Dyn(:,2) = Dyn(:,2) - Dyn(1,2);
Dyn(:,3) = Dyn(:,3) - Dyn(1,3);
Dyn(:,4) = Dyn(:,4) - Dyn(1,4);
Dyn(:,5) = Dyn(:,5) - Dyn(1,5);
Dyn(:,6) = Dyn(:,6) - Dyn(1,6);
Dyn = -Dyn;

frequency = zeros(size(Dyn, 1), 3);
frequency_fcc_cal = zeros(size(Dyn, 1), 1);
for k = 2:size(Dyn, 1)      % 第一个[0,0,0]不要考虑
    temp = [Dyn(k, 1), Dyn(k, 4), Dyn(k, 6);
        Dyn(k, 4), Dyn(k, 2), Dyn(k, 5);
        Dyn(k, 6), Dyn(k, 5), Dyn(k, 3)];
    frequency(k, :) = eigs(temp);
    frequency_fcc_cal(k) = 10*sqrt(max(frequency(k, :)))/2/pi;
end
end