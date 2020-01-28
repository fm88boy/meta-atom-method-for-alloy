%% 优化程序

function A_Fitting()
initialize();         % 设定随机数与子函数路径
[x0, para, lb, ub, A, b, options, MAX_iter, rand_ratio] = read_parameter();

%% 优化算法
for cons_iter = 1:MAX_iter  % 循环MAX_iter次
    disp('Start optimizing...')
    [x, fval] = fmincon(@(data)Inner_fun(data, para), x0, A, b, [], [], lb, ub, [], options);
    [out_txt, ~] = Post_analysis(x, para, 0, {});   % 根据版本号子目录下的后处理程序来分析

    Out_converge(out_txt, x, fval, cons_iter);
    name = sprintf('Output_%d', cons_iter);
    save(name, 'x', 'para')
    
    load('Output_lowest.mat', 'x', 'para')        % 下一次修改应该是从最好的数据来出发，而不应该是从比较差的数据出发。
    [x0, para] = random_x_and_weight(x, para, rand_ratio);
end
end

function initialize()
number_processors = 12;     % 调用的CPU数
sub_path = sprintf('%s/SubFunction', pwd);
addpath(sub_path);

fprintf('Date:\t%s\n\n', datestr(now, 31))      % 用时间的秒数来设定随机数
seed = str2double(datestr(now, 'FFFHHMMSS'));
s = RandStream('mt19937ar', 'Seed', seed);
RandStream.setGlobalStream(s)
fprintf('Random seed: %d\n', seed)

wait_time = 100*rand(1);
fprintf('Pause %f second\n', wait_time)

if ~ispc        % 如果不是Windows
    pause(wait_time)
    try
        parpool(number_processors);
    catch
        pause(30)
        parpool(number_processors);
    end
end
end

function [x, para, lb, ub, A, b, options, MAX_iter, rand_ratio] = read_parameter()
coeff = [1e-6, 10];   % 第1项为delta步长，第2项为调整步长
tol = [1e-6, 1e-3];         % 对于优化函数收敛的条件，第二个是对于cons的约束条件
MAX_iter = 500;         % 最多允许的约束循环次数
MAX_step = 800;        % 单次优化最大步数
MAX_cons = 10000;       % 最多的约束和目标量
rand_ratio = 0.8;   % 每次对原来的值扰动的最大比例

% 在后续的代码中，如果有三项相连，则依次为FCC、HCP和BCC相
weight.elastic = [5, 0, 0];         % 弹性常数
weight.dyn_lattice_energy = [1, 1, 1.5];     % 动态晶格常数与能量
weight.surf_energy = [1, 0, 0];     % 表面能
weight.sf_energy = [1, 0, 0];       % 层错能与不稳定层错能
weight.vac_energy = [0.5, 0, 0];    % 空穴能
weight.fcc_strain = [1, 1, 1];      % FCC应变路径能量
weight.prevent_hcp_energy = 0;      % 确保HCP相最低能量比FCC相高的能量
weight.phonon = [0, 0, 0];          % 声子
weight.phi_well = 1;                % 对于phi函数的势肼深度的优化
weight.fcc_specific_config = 1;     % 对于特定构型的限制（构型---能量）

weight = norm_struct(weight);       % 再将其归一化

% 单位说明：如果是长度，单位为Ang，如果为能量，单位为eV，如果界面能，单位为mJ/m2，如果是弹性常数，单位为GPa

cons_data.lattice = [3.604, [2.481, 3.933], 2.820];             % 晶格常数，依次为FCC、HCP_a、HCP_c和BCC
cons_data.atomic_weight = 55.4968;              % 平均的原子质量
cons_data.fcc_elastic = [417, 201, 329];        % FCC弹性常数，按照顺序，依次为C11,C12,C44
cons_data.hcp_elastic = [0, 0, 0, 0, 0];        % HCP弹性常数，按照顺序，依次为C11,C12,C13,C33,C44
cons_data.bcc_elastic = [0, 0, 0];              % BCC弹性常数，按照顺序，依次为C11,C12,C44
cons_data.sub_energy = [4.745, 4.788, 4.719];     % 升华能
cons_data.sf_energy = [-395, 0, 0];             % 层错能
cons_data.usf_energy = [186, 0, 0];             % 不稳定层错能
cons_data.fcc_surface_energy = [2367, 2718, 2783];      % FCC表面能，依次为111,110和100面
cons_data.hcp_surface_energy = [0, 0, 0];       % HCP表面能，依次为111,110和100面
cons_data.bcc_surface_energy = [0, 0, 0];       % BCC表面能，依次为111,110和100面
cons_data.vac_energy = [6.52, 0, 0];          % 空穴形成能
cons_data.fcc_q_mesh = [0.0, 0.0, 0.0, 0.0];    % FCC的Q点和目标频率THz
cons_data.hcp_q_mesh = [0.0, 0.0, 0.0, 0.0];    % HCP的Q点和目标频率THz
cons_data.bcc_q_mesh = [0.0, 0.0, 0.0, 0.0];    % BCC的Q点和目标频率THz
cons_data.fcc_strain_100 = [ -0.15, 0.05, 1.0;
                             -0.05, 0.03, 1.0;
                              0.05, 0.03, 1.0;
                              0.15, 0.05, 1.0];     % 100路径的应变（Bain路径），目标能量最小值（eV），目标能量最大值（eV）
cons_data.fcc_strain_110 = [ -0.05, 0.02, 1.0;
                              0.05, 0.02, 1.0;];    % 110路径的应变，目标能量最小值（eV），目标能量最大值（eV）
cons_data.fcc_strain_111 = [ -0.05, 0.02, 1.0;
                              0.05, 0.02, 1.0;];    % 111路径的应变，目标能量最小值（eV），目标能量最大值（eV）
cons_data.prevent_hcp_energy = [-0.05, 0.05, -1.0, 1.0];    % HCP能量扫描，z轴应变下区间与上区间，目标能量的下区间和上区间（eV）
cons_data.phi_well = [2.2, 4.0, 5, 0.4, 2.3, 10]; % phi函数形状的设置：分割点1，分割点2，分割点3，区间1与区间2的最小差距
                                                  % 最大值控制范围的下限，最大许可phi值

cons_data.fcc_elastic = cons_data.fcc_elastic * 1/160;      % 量纲转化
cons_data.hcp_elastic = cons_data.hcp_elastic * 1/160;
cons_data.bcc_elastic = cons_data.bcc_elastic * 1/160;
cons_data.sf_energy = cons_data.sf_energy * 1/16*1e-3;
cons_data.usf_energy = cons_data.usf_energy * 1/16*1e-3;
cons_data.fcc_surface_energy = cons_data.fcc_surface_energy * 1/16*1e-3;
cons_data.hcp_surface_energy = cons_data.hcp_surface_energy * 1/16*1e-3;
cons_data.bcc_surface_energy = cons_data.bcc_surface_energy * 1/16*1e-3;

control.F_phi_range = [500, 500];                   % F函数自变量的取值范围，phi自变量的取值范围，F(1)可调节范围
control.fcc_relation_C12_C44 = [1, 0.5];            % 是否限制C12-C44<0.5*(C12_cons-C44_cons)确保C12与C44的距离（仅适用于C12<C44），比例系数
control.dF_restrict = [1, 6, 5, 40];                % 是否限制F'<0，采样点数，生效范围下区间，生效范围上区间
control.fcc_dF_restrict_equal = [1, -0.1, 3];       % 是否限制dF_rho_fcc值，限制的最大值，附近的采样点数
control.phi_restrict = [1, 100];                    % 是否限制phi(knots(1))的值，限制的最大值
control.dphi_restrict = [1, 1.8, 2.0, 2.1];         % 是否限制phi'<0，施加限制的位置1，2，3
control.ddphi_restrict = [1, 3, 10, 2.5, 4.5];      % 是否限制phi''的绝对值，限制的绝对值，采样点数，生效范围下区间，生效范围上区间
control.fcc_rose_equation = [1, 11, 0.15, 1e-4, 2];    % 是否限制FCC的rose方程，采样点数，扫描范围的绝对值，底部限制，约束的模式
control.fcc_100_tension = [1, 0.05, 0.6, -0.06, 0.6, 0.12, 0.4]; 
% 是否限制FCC的100方向的拉伸，应变点1，最少能量百分比1，应变点2，最少能量百分比2
control.hcp_rose_equation = [1, 11, 0.15, 1e-4, 2];    % 是否限制HCP的rose方程，采样点数，a扫描范围绝对值，a底部限制，约束的模式
control.hcp_a_c_range = [1, 0.05, 0.6, -0.06, 0.6, 0.12, 0.4]; 
control.bcc_rose_equation = [1, 6, 0.15, 1e-4, 1];    % 是否限制BCC的rose方程，采样点数，a扫描范围绝对值，a底部限制，约束的模式
control.fcc_sf_relation = 1;                        % 是否限制FCC的GPFE一般形式，不稳定层错能为最大值（4个不等式），并且Esf<Eusf
control.fcc_surface_relation = 1;                   % 是否限制E_fcc_111<E_fcc_100<E_fcc_110

% cutoff = cons_data(1) * 1.65;    %比usf所需的数据略微大了一点点
cutoff = ceil(cons_data.lattice(1) * 1.65 * 10)/10;     % 根据FCC晶格选取cutoff半径
range = [0, 300, 0, cutoff]; % 第一个是F函数的取值范围，第二个是phi函数的取值范围

load start_point x para
F_data = x(1:para.split-1);   % 把数据，分为F函数与phi函数
phi_data = x(para.split:end);
para = replot_rho(para, range);      % 确保rho函数的取值范围正确，消除rho函数中的NaN

para.range = range;
para.cons_data = cons_data;
para.weight = weight;

if weight.fcc_specific_config ~= 0       % 如果需要读入特殊构型
    para = load_fcc_specific_config(para);
end

F_phi_cons = F_phi_cons_create(MAX_cons);    % 给出保存约束条件的结构体
[x, lb, ub, A, b, F_phi_cons] = make_cons(F_data, phi_data, para, cons_data, control, F_phi_cons);
F_phi_cons = F_phi_cons_target(para, F_phi_cons);        % 根据优化目标量，写入线性规划的模块
F_phi_cons = F_phi_cons_regulate(F_phi_cons);            % 将约束条件缩并且规整

A = A .* repmat(para.base_data, size(A, 1), 1);     % 将百分比乘进去
[para.DynData, para.cons_q] = select_fcc_q_mesh(cons_data.fcc_q_mesh, para.q_list, para.D_list);

fval = 100000;
cons_iter = 1;      % 用来标记这个计算是第几步
[x, para] = random_x_and_weight(x, para, rand_ratio);     % 权重也可以自动调整

save Output_lowest x para fval cons_iter

options = optimset('Algorithm', 'active-set', 'DiffMaxChange', coeff(2), 'DiffMinChange', coeff(1), 'Display', 'iter-detailed', ...
    'FinDiffType', 'forward', 'FunValCheck', 'on', 'MaxFunEvals', 1e8, 'MaxIter', MAX_step, 'OutputFcn', @Output, ...
    'TolFun', tol(1), 'TolX', tol(1), 'TolCon', tol(2), 'UseParallel', ~ispc, 'RelLineSrchBnd', coeff(2), 'RelLineSrchBndDuration', 100000);
% 如果是windows就不要并行了
end

function out = Inner_fun(data, para)
% 建立起F，rho和phi函数
rho = para.rho;
weight = para.weight;
critical = para.critical;   % 弹性常数所需要的几何因子
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

data = data .* para.base_data;      % 得到基础值
[F, phi] = F_phi(para.range, data, para.split, para.knots, para.order, 10000);

a_fcc = para.cons_data.lattice(1);
a_hcp = para.cons_data.lattice(2);
c_hcp = para.cons_data.lattice(3);
a_bcc = para.cons_data.lattice(4);
C11_fcc = para.cons_data.fcc_elastic(1);
C12_fcc = para.cons_data.fcc_elastic(2);
C44_fcc = para.cons_data.fcc_elastic(3);
Esub_fcc = para.cons_data.sub_energy(1);
Esub_hcp = para.cons_data.sub_energy(2);
Esub_bcc = para.cons_data.sub_energy(3);
esf_fcc = para.cons_data.sf_energy(1);
eusf_fcc = para.cons_data.usf_energy(1);
e111_fcc = para.cons_data.fcc_surface_energy(1);
e110_fcc = para.cons_data.fcc_surface_energy(2);
e100_fcc = para.cons_data.fcc_surface_energy(3);
evac_fcc = para.cons_data.vac_energy(1);
fcc_strain_100 = para.cons_data.fcc_strain_100;
fcc_strain_110 = para.cons_data.fcc_strain_110;
fcc_strain_111 = para.cons_data.fcc_strain_111;
fcc_specific_config = para.fcc_specific_config;
fcc_specific_energy = para.fcc_specific_energy;
prevent_hcp_energy = para.cons_data.prevent_hcp_energy;
phi_well = para.cons_data.phi_well;
cons_q = para.cons_q;

drho = (rho(3:end, 2) - rho(1:end-2, 2)) ./ (rho(3:end, 1) - rho(1:end-2, 1));
drho_x = rho(2:end-1, 1);
ddrho = (drho(3:end) - drho(1:end-2)) ./ (drho_x(3:end) - drho_x(1:end-2));
ddrho_x = drho_x(2:end-1);

dphi = (phi(3:end, 2) - phi(1:end-2, 2)) ./ (phi(3:end, 1) - phi(1:end-2, 1));
dphi_x = phi(2:end-1, 1);
ddphi = (dphi(3:end) - dphi(1:end-2)) ./ (dphi_x(3:end) - dphi_x(1:end-2));
ddphi_x = dphi_x(2:end-1);

dF = (F(3:end, 2) - F(1:end-2, 2)) ./ (F(3:end, 1) - F(1:end-2, 1));
dF_x = F(2:end-1, 1);
ddF = (dF(3:end) - dF(1:end-2)) ./ (dF_x(3:end) - dF_x(1:end-2));
ddF_x = dF_x(2:end-1);

out = 0;

%% 优化FCC的弹性常数
if weight.elastic(1) ~= 0
    r_fcc = alpha_fcc*a_fcc;
    rho_fcc = sum(N_fcc .* interp1(rho(:,1), rho(:,2), r_fcc));

    omega = a_fcc^3/4;      % 单个原子所占体积，这里面使用FCC的晶格常数来计算
    drho_fcc = interp1(drho_x, drho, r_fcc);        % 求出fcc相的各个导数
    ddrho_fcc = interp1(ddrho_x, ddrho, r_fcc);
    dphi_fcc = interp1(dphi_x, dphi, r_fcc);
    ddphi_fcc = interp1(ddphi_x, ddphi, r_fcc);
    dF_fcc = interp1(dF_x, dF, rho_fcc);
    ddF_fcc = interp1(ddF_x, ddF, rho_fcc);
    
    [C11_fcc_cal, C12_fcc_cal, C44_fcc_cal] = fit_fcc_elastic(ddF_fcc, dF_fcc, ddphi_fcc, dphi_fcc, ddrho_fcc, drho_fcc, ...
        alpha_fcc, a_fcc, critical, omega);
    
    out = out + weight.elastic(1)/3 * (C11_fcc - C11_fcc_cal)^2/C11_fcc^2;
    out = out + weight.elastic(1)/3 * (C12_fcc - C12_fcc_cal)^2/C12_fcc^2;
    out = out + weight.elastic(1)/3 * (C44_fcc - C44_fcc_cal)^2/C44_fcc^2;
end

%% 优化HCP的弹性常数，C11,C12,C13,C33,C44
if weight.elastic(2) ~= 0
    % 等待编写纯HCP的弹性常数方程
end

%% 优化BCC的弹性常数，C11,C12,C44
if weight.elastic(3) ~= 0
    % 等待编写纯BCC的弹性常数方程
end

%% 优化FCC、HCP和BCC的晶格常数与内聚能
if weight.dyn_lattice_energy(1) ~= 0    % 动态测量，FCC的a和能量
    weight_single = 0.5 * weight.dyn_lattice_energy(1);
    [a_fcc_cal, fcc_energy_cal] = fit_fcc_real_a_and_energy(F, rho, phi, a_fcc, para.range(4));
    out = out + weight_single * (a_fcc_cal - a_fcc)^2/(a_fcc)^2;
    out = out + weight_single * (abs(fcc_energy_cal) - abs(Esub_fcc))^2/(Esub_fcc)^2;
end

if weight.dyn_lattice_energy(2) ~= 0    % 动态测量，HCP的a,c和能量
    weight_single = 0.5 * weight.dyn_lattice_energy(2);
    [a_hcp_cal, c_hcp_cal, hcp_energy_cal] = fit_hcp_real_a_c_and_energy(F, rho, phi, a_hcp, c_hcp, para.range(4));
    out = out + 0.5 * weight_single * (a_hcp_cal - a_hcp)^2/(a_hcp)^2;
    out = out + 0.5 * weight_single * (c_hcp_cal - c_hcp)^2/(c_hcp)^2;
    out = out + weight_single * (abs(hcp_energy_cal) - abs(Esub_hcp))^2/(Esub_hcp)^2;
end

if weight.dyn_lattice_energy(3) ~= 0    % 动态测量，BCC的a和能量
    weight_single = 0.5 * weight.dyn_lattice_energy(3);
    [a_bcc_cal, bcc_energy_cal] = fit_bcc_real_a_and_energy(F, rho, phi, a_bcc, para.range(4));
    out = out + weight_single * (a_bcc_cal - a_bcc)^2/(a_bcc)^2;
    out = out + weight_single * (abs(bcc_energy_cal) - abs(Esub_bcc))^2/(Esub_bcc)^2;
end

%% 优化FCC的能量
% FCC的层错能
if weight.sf_energy(1) ~= 0     % FCC的层错能
    [esf_fcc_cal, eusf_fcc_cal] = fit_fcc_sf(F, rho, phi, alpha_fcc_sf, N_fcc_sf, alpha_fcc_usf1, N_fcc_usf1, ...
        alpha_fcc_usf2, N_fcc_usf2, a_fcc, -1*Esub_fcc);
    out = out + weight.sf_energy(1) * 0.7 * (esf_fcc - esf_fcc_cal)^2/esf_fcc^2;        % 对于层错能，给了60%的权重
    out = out + weight.sf_energy(1) * 0.3 * (eusf_fcc - eusf_fcc_cal)^2/eusf_fcc^2;     % 对于不稳定层错能，给了40%的权重
end

% FCC的111,110,100表面能
if weight.surf_energy(1) ~= 0
    [e111_fcc_cal, e110_fcc_cal, e100_fcc_cal] = fit_fcc_surface(F, rho, phi, N_fcc_111, N_fcc_110, N_fcc_100, r_fcc, a_fcc, -1*Esub_fcc);
    out = out + weight.surf_energy(1)/3 * (e111_fcc - e111_fcc_cal)^2/e111_fcc^2;
    out = out + weight.surf_energy(1)/3 * (e110_fcc - e110_fcc_cal)^2/e110_fcc^2;
    out = out + weight.surf_energy(1)/3 * (e100_fcc - e100_fcc_cal)^2/e100_fcc^2;
end

% FCC的空穴能
if weight.vac_energy(1) ~= 0
    evac_fcc_cal = fit_fcc_vac(F, rho, phi, N_fcc_vac, N_fcc, -1*Esub_fcc, r_fcc);    
    out = out + weight.vac_energy(1) * (evac_fcc - evac_fcc_cal)^2/evac_fcc^2;
end

% FCC沿100方向的拉伸，Bain路径
if weight.fcc_strain(1) ~= 0
    fcc_strain_100_cal = fit_fcc_strain_100(F, rho, phi, fcc_strain_100, a_fcc, Esub_fcc, para.range(4));
    weight_single = weight.fcc_strain(1)/size(fcc_strain_100_cal, 1);
    
    for i = 1:size(fcc_strain_100_cal, 1)
        if fcc_strain_100_cal(i, 2) < fcc_strain_100(i, 2)    % 如果能量并没有落入区间内，则给予惩罚
            out = out + weight_single * (fcc_strain_100(i, 2) - fcc_strain_100_cal(i, 2))^2/fcc_strain_100(i, 2)^2;
        elseif fcc_strain_100_cal(i, 2) > fcc_strain_100(i, 3)
            out = out + weight_single * (fcc_strain_100(i, 3) - fcc_strain_100_cal(i, 2))^2/fcc_strain_100(i, 3)^2;
        end
    end
end

% FCC沿110方向的拉伸
if weight.fcc_strain(2) ~= 0
    fcc_strain_110_cal = fit_fcc_strain_110(F, rho, phi, fcc_strain_110, a_fcc, Esub_fcc, para.range(4));
    weight_single = weight.fcc_strain(2)/size(fcc_strain_110_cal, 1);
    
    for i = 1:size(fcc_strain_110_cal, 1)
        if fcc_strain_110_cal(i, 2) < fcc_strain_110(i, 2)    % 如果能量并没有落入区间内，则给予惩罚
            out = out + weight_single * (fcc_strain_110(i, 2) - fcc_strain_110_cal(i, 2))^2/fcc_strain_110(i, 2)^2;
        elseif fcc_strain_110_cal(i, 2) > fcc_strain_110(i, 3)
            out = out + weight_single * (fcc_strain_110(i, 3) - fcc_strain_110_cal(i, 2))^2/fcc_strain_110(i, 3)^2;
        end
    end
end

% FCC沿111方向的拉伸
if weight.fcc_strain(3) ~= 0
    fcc_strain_111_cal = fit_fcc_strain_111(F, rho, phi, fcc_strain_111, a_fcc, Esub_fcc, para.range(4));
    weight_single = weight.fcc_strain(3)/size(fcc_strain_111_cal, 1);
    
    for i = 1:size(fcc_strain_111_cal, 1)
        if fcc_strain_111_cal(i, 2) < fcc_strain_111(i, 2)    % 如果能量并没有落入区间内，则给予惩罚
            out = out + weight_single * (fcc_strain_111(i, 2) - fcc_strain_111_cal(i, 2))^2/fcc_strain_111(i, 2)^2;
        elseif fcc_strain_111_cal(i, 2) > fcc_strain_111(i, 3)
            out = out + weight_single * (fcc_strain_111(i, 3) - fcc_strain_111_cal(i, 2))^2/fcc_strain_111(i, 3)^2;
        end
    end
end

%% 对于特定构型的限制（构型---能量）
if weight.fcc_specific_config ~= 0
    weight_single = weight.fcc_specific_config/size(fcc_specific_config.geo, 1);
    fcc_specific_energy_cal = fit_fcc_specific_config(F, rho, phi, a_fcc, fcc_specific_config);
    
    for i = 1:size(fcc_specific_energy_cal, 1)
        if fcc_specific_energy_cal(i, 1) < fcc_specific_energy(i, 1)
            out = out + weight_single * (fcc_specific_energy_cal(i, 1) - fcc_specific_energy(i, 1))^2/fcc_specific_energy(i, 1)^2;
        elseif fcc_specific_energy_cal(i, 1) > fcc_specific_energy(i, 2)
            out = out + weight_single * (fcc_specific_energy_cal(i, 1) - fcc_specific_energy(i, 2))^2/fcc_specific_energy(i, 2)^2;
        end
    end
end

%% HCP与FCC的关系
if weight.prevent_hcp_energy ~= 0
    hcp_energy = fit_prevent_hcp(F, rho, phi, prevent_hcp_energy(1:2), a_fcc, Esub_fcc, para.range(4));
    if hcp_energy < prevent_hcp_energy(3)
        out = out + weight.prevent_hcp_energy * (hcp_energy - prevent_hcp_energy(3))^2/0.01^2;  % 这里面选一个较大的量用来判断
    elseif hcp_energy > prevent_hcp_energy(4)
        out = out + weight.prevent_hcp_energy * (hcp_energy - prevent_hcp_energy(4))^2/prevent_hcp_energy(4)^2;
    end
end

%% phi函数的形状
if weight.phi_well ~= 0
    [diff_cal, max_phi_cal] = fit_phi_well(phi, phi_well);
    if diff_cal < phi_well(4)       % 这个值控制了第一个势肼相对第二个势肼的深度
        out = out + weight.phi_well/2 * abs(diff_cal - phi_well(4));
    end
    if max_phi_cal > phi_well(6)
        out = out + weight.phi_well/2 * abs(max_phi_cal - phi_well(6));
    end
end

%% FCC中Q点的频率
if weight.phonon(1) ~= 0
    frequency_fcc_cal = fit_fcc_frequency(para.DynData, ddphi_fcc, dphi_fcc, ddF_fcc, dF_fcc, ddrho_fcc, drho_fcc, alpha_fcc, r_fcc);
    weight_single = weight.phonon(1)/(size(cons_q, 1)-1);  % 计算出单个变量的权重
    out = out + sum(weight_single * (frequency_fcc_cal - cons_q).^2./cons_q.^2);
end
end

function stop = Output(x, optimValues, state)
stop = false;
load Output_lowest para
fval = optimValues.fval;
switch state
    case 'iter'
        switch mod(optimValues.iteration, 2000)      % 每1000步做一次IO操作，这样可以加速计算速度
            case 0
                save Output_iter_1 x para fval;
            case 1000
                save Output_iter_2 x para fval;
        end
end
end
