%% 根据优化的目标量，写入线性规划的条件

function F_phi_cons = F_phi_cons_target(para, F_phi_cons)
% 建立起F，rho和phi函数
rho = para.rho;         % 记录了电荷密度函数
r_cut = para.range(4);
weight = para.weight;
alpha_fcc = para.alpha_fcc;
N_fcc = para.N_fcc;
alpha_fcc_sf = para.alpha_fcc_sf;
N_fcc_sf = para.N_fcc_sf;
alpha_fcc_usf1 = para.alpha_fcc_usf1;
N_fcc_usf1 = para.N_fcc_usf1;
alpha_fcc_usf2 = para.alpha_fcc_usf2;
N_fcc_usf2 = para.N_fcc_usf2;
N_fcc_111 = para.N_fcc_111;
N_fcc_110 = para.N_fcc_110;
N_fcc_100 = para.N_fcc_100;
N_fcc_vac = para.N_fcc_vac;

a_fcc = para.cons_data.lattice(1);
a_hcp = para.cons_data.lattice(2);
c_hcp = para.cons_data.lattice(3);
a_bcc = para.cons_data.lattice(4);
Esub_fcc = para.cons_data.sub_energy(1);
Esub_hcp = para.cons_data.sub_energy(2);
Esub_bcc = para.cons_data.sub_energy(3);
esf_fcc = para.cons_data.sf_energy(1);
eusf_fcc = para.cons_data.usf_energy(1);
e111_fcc = para.cons_data.fcc_surface_energy(1);
e110_fcc = para.cons_data.fcc_surface_energy(2);
e100_fcc = para.cons_data.fcc_surface_energy(3);
evac_fcc = para.cons_data.vac_energy(1);

%% 优化FCC、HCP和BCC的晶格常数与内聚能
if weight.dyn_lattice_energy(1) ~= 0               % FCC的a和能量
    radius_fcc = alpha_fcc * a_fcc;
    rho_fcc = sum(N_fcc .* interp1(rho(:,1), rho(:,2), radius_fcc));
    F_coeff = [1, rho_fcc];
    phi_coeff = [0.5*N_fcc', radius_fcc'];
    b_coeff = [2, -Esub_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
end

if weight.dyn_lattice_energy(2) ~= 0                % HCP的a,c和能量
    [radius_hcp, N_hcp] = cal_hcp_radius(a_hcp, c_hcp, r_cut);         % 获得新的HCP半径和原子数
    rho_hcp = sum(N_hcp .* interp1(rho(:,1), rho(:,2), radius_hcp));
    F_coeff = [1, rho_hcp];
    phi_coeff = [0.5*N_hcp', radius_hcp'];
    b_coeff = [2, -Esub_hcp];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
end

if weight.dyn_lattice_energy(3) ~= 0    % 动态测量，BCC的a和能量
    [radius_bcc, N_bcc] = cal_bcc_radius(a_bcc, r_cut);      % 计算出对应的几何因子
    rho_bcc = sum(N_bcc .* interp1(rho(:,1), rho(:,2), radius_bcc));   % 计算出对应的电荷密度
    F_coeff = [1, rho_bcc];
    phi_coeff = [0.5*N_bcc', radius_bcc'];
    b_coeff = [2, -Esub_bcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
end

%% 优化FCC的能量
% FCC的层错能
if weight.sf_energy(1) ~= 0     % FCC的层错能
    radius_fcc_sf = alpha_fcc_sf * a_fcc;
    rho_fcc_sf = sum(N_fcc_sf .* interp1(rho(:,1), rho(:,2), radius_fcc_sf));
    F_coeff = [1, rho_fcc_sf;
               -1, rho_fcc];         % 必须要减去完美的FCC相的数据
    phi_coeff = [0.5*N_fcc_sf', radius_fcc_sf';
                 -0.5*N_fcc', radius_fcc'];
    F_coeff(:, 1) = 16*sqrt(3)/3/a_fcc^2 * F_coeff(:, 1);
    phi_coeff(:, 1) = 16*sqrt(3)/3/a_fcc^2 * phi_coeff(:, 1);
    b_coeff = [2, esf_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);

    radius_fcc_usf1 = alpha_fcc_usf1 * a_fcc;
    radius_fcc_usf2 = alpha_fcc_usf2 * a_fcc;
    rho_fcc_usf1 = sum(N_fcc_usf1 .* interp1(rho(:,1), rho(:,2), radius_fcc_usf1));
    rho_fcc_usf2 = sum(N_fcc_usf2 .* interp1(rho(:,1), rho(:,2), radius_fcc_usf2));
    F_coeff = [1, rho_fcc_usf1; 
               1, rho_fcc_usf2;
               -2, rho_fcc];         % 必须要减去完美的FCC相的数据
    phi_coeff = [0.5*N_fcc_usf1', radius_fcc_usf1';
                 0.5*N_fcc_usf2', radius_fcc_usf2';
                 -2*0.5*N_fcc', radius_fcc'];
    F_coeff(:, 1) = 8*sqrt(3)/3/a_fcc^2 * F_coeff(:, 1);
    phi_coeff(:, 1) = 8*sqrt(3)/3/a_fcc^2 * phi_coeff(:, 1);
    b_coeff = [2, eusf_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
end

% FCC的111,110,100表面能
if weight.surf_energy(1) ~= 0
    rho_fcc_111_1 = sum(N_fcc_111(1, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_111_2 = sum(N_fcc_111(2, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    F_coeff = [1, rho_fcc_111_1;
               1, rho_fcc_111_2;
               -2, rho_fcc];         % 必须要减去完美的FCC相的数据
    phi_coeff = [0.5*N_fcc_111(1, :)', radius_fcc';
                 0.5*N_fcc_111(2, :)', radius_fcc';
                 -2*0.5*N_fcc', radius_fcc'];
    F_coeff(:, 1) = 4*sqrt(3)/3/a_fcc^2 * F_coeff(:, 1);
    phi_coeff(:, 1) = 4*sqrt(3)/3/a_fcc^2 * phi_coeff(:, 1);
    b_coeff = [2, e111_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    
    rho_fcc_110_1 = sum(N_fcc_110(1, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_110_2 = sum(N_fcc_110(2, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_110_3 = sum(N_fcc_110(3, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_110_4 = sum(N_fcc_110(4, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    F_coeff = [1, rho_fcc_110_1;
               1, rho_fcc_110_2;
               1, rho_fcc_110_3;
               1, rho_fcc_110_4;
               -4, rho_fcc];         % 必须要减去完美的FCC相的数据
    phi_coeff = [0.5*N_fcc_110(1, :)', radius_fcc';
                 0.5*N_fcc_110(2, :)', radius_fcc';
                 0.5*N_fcc_110(3, :)', radius_fcc';
                 0.5*N_fcc_110(4, :)', radius_fcc';
                 -4*0.5*N_fcc', radius_fcc'];
    F_coeff(:, 1) = sqrt(2)/a_fcc^2 * F_coeff(:, 1);
    phi_coeff(:, 1) = sqrt(2)/a_fcc^2 * phi_coeff(:, 1);
    b_coeff = [2, e110_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);

    rho_fcc_100_1 = sum(N_fcc_100(1, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_100_2 = sum(N_fcc_100(2, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_100_3 = sum(N_fcc_100(3, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    F_coeff = [1, rho_fcc_100_1;
               1, rho_fcc_100_2;
               1, rho_fcc_100_3;
               -3, rho_fcc];         % 必须要减去完美的FCC相的数据
    phi_coeff = [0.5*N_fcc_100(1, :)', radius_fcc';
                 0.5*N_fcc_100(2, :)', radius_fcc';
                 0.5*N_fcc_100(3, :)', radius_fcc';
                 -3*0.5*N_fcc', radius_fcc'];
    F_coeff(:, 1) = 2/a_fcc^2 * F_coeff(:, 1);
    phi_coeff(:, 1) = 2/a_fcc^2 * phi_coeff(:, 1);
    b_coeff = [2, e100_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
end

% FCC的空穴能
if weight.vac_energy(1) ~= 0
    rho_fcc_vac_1 = sum(N_fcc_vac(1, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_vac_2 = sum(N_fcc_vac(2, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_vac_3 = sum(N_fcc_vac(3, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_vac_4 = sum(N_fcc_vac(4, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    rho_fcc_vac_5 = sum(N_fcc_vac(5, :) .* interp1(rho(:,1), rho(:,2), radius_fcc));
    F_coeff = [N_fcc(1), rho_fcc_vac_1;
               N_fcc(2), rho_fcc_vac_2;
               N_fcc(3), rho_fcc_vac_3;
               N_fcc(4), rho_fcc_vac_4;
               N_fcc(5), rho_fcc_vac_5;
               -sum(N_fcc), rho_fcc];         % 必须要减去完美的FCC相的数据
    phi_coeff = [N_fcc(1)*0.5*N_fcc_vac(1, :)', radius_fcc';
                 N_fcc(2)*0.5*N_fcc_vac(2, :)', radius_fcc';
                 N_fcc(3)*0.5*N_fcc_vac(3, :)', radius_fcc';
                 N_fcc(4)*0.5*N_fcc_vac(4, :)', radius_fcc';
                 N_fcc(5)*0.5*N_fcc_vac(5, :)', radius_fcc';
                 -sum(N_fcc)*0.5*N_fcc', radius_fcc'];
    b_coeff = [2, evac_fcc];
    F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
end

end