%% 根据FCC的几何因子计算系统能量

function energy = cal_fcc_geo_energy(F, rho, phi, a_fcc, geo, F_geo)
radius = geo(:, 1) * a_fcc;
phi_sum = sum(geo(:, 2) .* interp1(phi(:,1), phi(:,2), radius));

F_sum = 0;
for i = 1:size(F_geo.factor, 1)
    radius = F_geo.factor{i}(:, 1) * a_fcc;
    rho_this = sum(F_geo.factor{i}(:, 2) .* interp1(rho(:,1), rho(:,2), radius));
    F_sum = F_sum + F_geo.count(i) * interp1(F(:,1), F(:,2), rho_this);
end
energy = F_sum + 0.5*phi_sum;
end
