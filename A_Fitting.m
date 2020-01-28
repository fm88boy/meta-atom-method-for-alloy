%% �Ż�����

function A_Fitting()
initialize();         % �趨��������Ӻ���·��
[x0, para, lb, ub, A, b, options, MAX_iter, rand_ratio] = read_parameter();

%% �Ż��㷨
for cons_iter = 1:MAX_iter  % ѭ��MAX_iter��
    disp('Start optimizing...')
    [x, fval] = fmincon(@(data)Inner_fun(data, para), x0, A, b, [], [], lb, ub, [], options);
    [out_txt, ~] = Post_analysis(x, para, 0, {});   % ���ݰ汾����Ŀ¼�µĺ������������

    Out_converge(out_txt, x, fval, cons_iter);
    name = sprintf('Output_%d', cons_iter);
    save(name, 'x', 'para')
    
    load('Output_lowest.mat', 'x', 'para')        % ��һ���޸�Ӧ���Ǵ���õ�����������������Ӧ���ǴӱȽϲ�����ݳ�����
    [x0, para] = random_x_and_weight(x, para, rand_ratio);
end
end

function initialize()
number_processors = 12;     % ���õ�CPU��
sub_path = sprintf('%s/SubFunction', pwd);
addpath(sub_path);

fprintf('Date:\t%s\n\n', datestr(now, 31))      % ��ʱ����������趨�����
seed = str2double(datestr(now, 'FFFHHMMSS'));
s = RandStream('mt19937ar', 'Seed', seed);
RandStream.setGlobalStream(s)
fprintf('Random seed: %d\n', seed)

wait_time = 100*rand(1);
fprintf('Pause %f second\n', wait_time)

if ~ispc        % �������Windows
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
coeff = [1e-6, 10];   % ��1��Ϊdelta��������2��Ϊ��������
tol = [1e-6, 1e-3];         % �����Ż������������������ڶ����Ƕ���cons��Լ������
MAX_iter = 500;         % ��������Լ��ѭ������
MAX_step = 800;        % �����Ż������
MAX_cons = 10000;       % ����Լ����Ŀ����
rand_ratio = 0.8;   % ÿ�ζ�ԭ����ֵ�Ŷ���������

% �ں����Ĵ����У����������������������ΪFCC��HCP��BCC��
weight.elastic = [5, 0, 0];         % ���Գ���
weight.dyn_lattice_energy = [1, 1, 1.5];     % ��̬������������
weight.surf_energy = [1, 0, 0];     % ������
weight.sf_energy = [1, 0, 0];       % ������벻�ȶ������
weight.vac_energy = [0.5, 0, 0];    % ��Ѩ��
weight.fcc_strain = [1, 1, 1];      % FCCӦ��·������
weight.prevent_hcp_energy = 0;      % ȷ��HCP�����������FCC��ߵ�����
weight.phonon = [0, 0, 0];          % ����
weight.phi_well = 1;                % ����phi������������ȵ��Ż�
weight.fcc_specific_config = 1;     % �����ض����͵����ƣ�����---������

weight = norm_struct(weight);       % �ٽ����һ��

% ��λ˵��������ǳ��ȣ���λΪAng�����Ϊ��������λΪeV����������ܣ���λΪmJ/m2������ǵ��Գ�������λΪGPa

cons_data.lattice = [3.604, [2.481, 3.933], 2.820];             % ������������ΪFCC��HCP_a��HCP_c��BCC
cons_data.atomic_weight = 55.4968;              % ƽ����ԭ������
cons_data.fcc_elastic = [417, 201, 329];        % FCC���Գ���������˳������ΪC11,C12,C44
cons_data.hcp_elastic = [0, 0, 0, 0, 0];        % HCP���Գ���������˳������ΪC11,C12,C13,C33,C44
cons_data.bcc_elastic = [0, 0, 0];              % BCC���Գ���������˳������ΪC11,C12,C44
cons_data.sub_energy = [4.745, 4.788, 4.719];     % ������
cons_data.sf_energy = [-395, 0, 0];             % �����
cons_data.usf_energy = [186, 0, 0];             % ���ȶ������
cons_data.fcc_surface_energy = [2367, 2718, 2783];      % FCC�����ܣ�����Ϊ111,110��100��
cons_data.hcp_surface_energy = [0, 0, 0];       % HCP�����ܣ�����Ϊ111,110��100��
cons_data.bcc_surface_energy = [0, 0, 0];       % BCC�����ܣ�����Ϊ111,110��100��
cons_data.vac_energy = [6.52, 0, 0];          % ��Ѩ�γ���
cons_data.fcc_q_mesh = [0.0, 0.0, 0.0, 0.0];    % FCC��Q���Ŀ��Ƶ��THz
cons_data.hcp_q_mesh = [0.0, 0.0, 0.0, 0.0];    % HCP��Q���Ŀ��Ƶ��THz
cons_data.bcc_q_mesh = [0.0, 0.0, 0.0, 0.0];    % BCC��Q���Ŀ��Ƶ��THz
cons_data.fcc_strain_100 = [ -0.15, 0.05, 1.0;
                             -0.05, 0.03, 1.0;
                              0.05, 0.03, 1.0;
                              0.15, 0.05, 1.0];     % 100·����Ӧ�䣨Bain·������Ŀ��������Сֵ��eV����Ŀ���������ֵ��eV��
cons_data.fcc_strain_110 = [ -0.05, 0.02, 1.0;
                              0.05, 0.02, 1.0;];    % 110·����Ӧ�䣬Ŀ��������Сֵ��eV����Ŀ���������ֵ��eV��
cons_data.fcc_strain_111 = [ -0.05, 0.02, 1.0;
                              0.05, 0.02, 1.0;];    % 111·����Ӧ�䣬Ŀ��������Сֵ��eV����Ŀ���������ֵ��eV��
cons_data.prevent_hcp_energy = [-0.05, 0.05, -1.0, 1.0];    % HCP����ɨ�裬z��Ӧ���������������䣬Ŀ��������������������䣨eV��
cons_data.phi_well = [2.2, 4.0, 5, 0.4, 2.3, 10]; % phi������״�����ã��ָ��1���ָ��2���ָ��3������1������2����С���
                                                  % ���ֵ���Ʒ�Χ�����ޣ�������phiֵ

cons_data.fcc_elastic = cons_data.fcc_elastic * 1/160;      % ����ת��
cons_data.hcp_elastic = cons_data.hcp_elastic * 1/160;
cons_data.bcc_elastic = cons_data.bcc_elastic * 1/160;
cons_data.sf_energy = cons_data.sf_energy * 1/16*1e-3;
cons_data.usf_energy = cons_data.usf_energy * 1/16*1e-3;
cons_data.fcc_surface_energy = cons_data.fcc_surface_energy * 1/16*1e-3;
cons_data.hcp_surface_energy = cons_data.hcp_surface_energy * 1/16*1e-3;
cons_data.bcc_surface_energy = cons_data.bcc_surface_energy * 1/16*1e-3;

control.F_phi_range = [500, 500];                   % F�����Ա�����ȡֵ��Χ��phi�Ա�����ȡֵ��Χ��F(1)�ɵ��ڷ�Χ
control.fcc_relation_C12_C44 = [1, 0.5];            % �Ƿ�����C12-C44<0.5*(C12_cons-C44_cons)ȷ��C12��C44�ľ��루��������C12<C44��������ϵ��
control.dF_restrict = [1, 6, 5, 40];                % �Ƿ�����F'<0��������������Ч��Χ�����䣬��Ч��Χ������
control.fcc_dF_restrict_equal = [1, -0.1, 3];       % �Ƿ�����dF_rho_fccֵ�����Ƶ����ֵ�������Ĳ�������
control.phi_restrict = [1, 100];                    % �Ƿ�����phi(knots(1))��ֵ�����Ƶ����ֵ
control.dphi_restrict = [1, 1.8, 2.0, 2.1];         % �Ƿ�����phi'<0��ʩ�����Ƶ�λ��1��2��3
control.ddphi_restrict = [1, 3, 10, 2.5, 4.5];      % �Ƿ�����phi''�ľ���ֵ�����Ƶľ���ֵ��������������Ч��Χ�����䣬��Ч��Χ������
control.fcc_rose_equation = [1, 11, 0.15, 1e-4, 2];    % �Ƿ�����FCC��rose���̣�����������ɨ�跶Χ�ľ���ֵ���ײ����ƣ�Լ����ģʽ
control.fcc_100_tension = [1, 0.05, 0.6, -0.06, 0.6, 0.12, 0.4]; 
% �Ƿ�����FCC��100��������죬Ӧ���1�����������ٷֱ�1��Ӧ���2�����������ٷֱ�2
control.hcp_rose_equation = [1, 11, 0.15, 1e-4, 2];    % �Ƿ�����HCP��rose���̣�����������aɨ�跶Χ����ֵ��a�ײ����ƣ�Լ����ģʽ
control.hcp_a_c_range = [1, 0.05, 0.6, -0.06, 0.6, 0.12, 0.4]; 
control.bcc_rose_equation = [1, 6, 0.15, 1e-4, 1];    % �Ƿ�����BCC��rose���̣�����������aɨ�跶Χ����ֵ��a�ײ����ƣ�Լ����ģʽ
control.fcc_sf_relation = 1;                        % �Ƿ�����FCC��GPFEһ����ʽ�����ȶ������Ϊ���ֵ��4������ʽ��������Esf<Eusf
control.fcc_surface_relation = 1;                   % �Ƿ�����E_fcc_111<E_fcc_100<E_fcc_110

% cutoff = cons_data(1) * 1.65;    %��usf�����������΢����һ���
cutoff = ceil(cons_data.lattice(1) * 1.65 * 10)/10;     % ����FCC����ѡȡcutoff�뾶
range = [0, 300, 0, cutoff]; % ��һ����F������ȡֵ��Χ���ڶ�����phi������ȡֵ��Χ

load start_point x para
F_data = x(1:para.split-1);   % �����ݣ���ΪF������phi����
phi_data = x(para.split:end);
para = replot_rho(para, range);      % ȷ��rho������ȡֵ��Χ��ȷ������rho�����е�NaN

para.range = range;
para.cons_data = cons_data;
para.weight = weight;

if weight.fcc_specific_config ~= 0       % �����Ҫ�������⹹��
    para = load_fcc_specific_config(para);
end

F_phi_cons = F_phi_cons_create(MAX_cons);    % ��������Լ�������Ľṹ��
[x, lb, ub, A, b, F_phi_cons] = make_cons(F_data, phi_data, para, cons_data, control, F_phi_cons);
F_phi_cons = F_phi_cons_target(para, F_phi_cons);        % �����Ż�Ŀ������д�����Թ滮��ģ��
F_phi_cons = F_phi_cons_regulate(F_phi_cons);            % ��Լ�����������ҹ���

A = A .* repmat(para.base_data, size(A, 1), 1);     % ���ٷֱȳ˽�ȥ
[para.DynData, para.cons_q] = select_fcc_q_mesh(cons_data.fcc_q_mesh, para.q_list, para.D_list);

fval = 100000;
cons_iter = 1;      % ���������������ǵڼ���
[x, para] = random_x_and_weight(x, para, rand_ratio);     % Ȩ��Ҳ�����Զ�����

save Output_lowest x para fval cons_iter

options = optimset('Algorithm', 'active-set', 'DiffMaxChange', coeff(2), 'DiffMinChange', coeff(1), 'Display', 'iter-detailed', ...
    'FinDiffType', 'forward', 'FunValCheck', 'on', 'MaxFunEvals', 1e8, 'MaxIter', MAX_step, 'OutputFcn', @Output, ...
    'TolFun', tol(1), 'TolX', tol(1), 'TolCon', tol(2), 'UseParallel', ~ispc, 'RelLineSrchBnd', coeff(2), 'RelLineSrchBndDuration', 100000);
% �����windows�Ͳ�Ҫ������
end

function out = Inner_fun(data, para)
% ������F��rho��phi����
rho = para.rho;
weight = para.weight;
critical = para.critical;   % ���Գ�������Ҫ�ļ�������
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

data = data .* para.base_data;      % �õ�����ֵ
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

%% �Ż�FCC�ĵ��Գ���
if weight.elastic(1) ~= 0
    r_fcc = alpha_fcc*a_fcc;
    rho_fcc = sum(N_fcc .* interp1(rho(:,1), rho(:,2), r_fcc));

    omega = a_fcc^3/4;      % ����ԭ����ռ�����������ʹ��FCC�ľ�����������
    drho_fcc = interp1(drho_x, drho, r_fcc);        % ���fcc��ĸ�������
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

%% �Ż�HCP�ĵ��Գ�����C11,C12,C13,C33,C44
if weight.elastic(2) ~= 0
    % �ȴ���д��HCP�ĵ��Գ�������
end

%% �Ż�BCC�ĵ��Գ�����C11,C12,C44
if weight.elastic(3) ~= 0
    % �ȴ���д��BCC�ĵ��Գ�������
end

%% �Ż�FCC��HCP��BCC�ľ��������ھ���
if weight.dyn_lattice_energy(1) ~= 0    % ��̬������FCC��a������
    weight_single = 0.5 * weight.dyn_lattice_energy(1);
    [a_fcc_cal, fcc_energy_cal] = fit_fcc_real_a_and_energy(F, rho, phi, a_fcc, para.range(4));
    out = out + weight_single * (a_fcc_cal - a_fcc)^2/(a_fcc)^2;
    out = out + weight_single * (abs(fcc_energy_cal) - abs(Esub_fcc))^2/(Esub_fcc)^2;
end

if weight.dyn_lattice_energy(2) ~= 0    % ��̬������HCP��a,c������
    weight_single = 0.5 * weight.dyn_lattice_energy(2);
    [a_hcp_cal, c_hcp_cal, hcp_energy_cal] = fit_hcp_real_a_c_and_energy(F, rho, phi, a_hcp, c_hcp, para.range(4));
    out = out + 0.5 * weight_single * (a_hcp_cal - a_hcp)^2/(a_hcp)^2;
    out = out + 0.5 * weight_single * (c_hcp_cal - c_hcp)^2/(c_hcp)^2;
    out = out + weight_single * (abs(hcp_energy_cal) - abs(Esub_hcp))^2/(Esub_hcp)^2;
end

if weight.dyn_lattice_energy(3) ~= 0    % ��̬������BCC��a������
    weight_single = 0.5 * weight.dyn_lattice_energy(3);
    [a_bcc_cal, bcc_energy_cal] = fit_bcc_real_a_and_energy(F, rho, phi, a_bcc, para.range(4));
    out = out + weight_single * (a_bcc_cal - a_bcc)^2/(a_bcc)^2;
    out = out + weight_single * (abs(bcc_energy_cal) - abs(Esub_bcc))^2/(Esub_bcc)^2;
end

%% �Ż�FCC������
% FCC�Ĳ����
if weight.sf_energy(1) ~= 0     % FCC�Ĳ����
    [esf_fcc_cal, eusf_fcc_cal] = fit_fcc_sf(F, rho, phi, alpha_fcc_sf, N_fcc_sf, alpha_fcc_usf1, N_fcc_usf1, ...
        alpha_fcc_usf2, N_fcc_usf2, a_fcc, -1*Esub_fcc);
    out = out + weight.sf_energy(1) * 0.7 * (esf_fcc - esf_fcc_cal)^2/esf_fcc^2;        % ���ڲ���ܣ�����60%��Ȩ��
    out = out + weight.sf_energy(1) * 0.3 * (eusf_fcc - eusf_fcc_cal)^2/eusf_fcc^2;     % ���ڲ��ȶ�����ܣ�����40%��Ȩ��
end

% FCC��111,110,100������
if weight.surf_energy(1) ~= 0
    [e111_fcc_cal, e110_fcc_cal, e100_fcc_cal] = fit_fcc_surface(F, rho, phi, N_fcc_111, N_fcc_110, N_fcc_100, r_fcc, a_fcc, -1*Esub_fcc);
    out = out + weight.surf_energy(1)/3 * (e111_fcc - e111_fcc_cal)^2/e111_fcc^2;
    out = out + weight.surf_energy(1)/3 * (e110_fcc - e110_fcc_cal)^2/e110_fcc^2;
    out = out + weight.surf_energy(1)/3 * (e100_fcc - e100_fcc_cal)^2/e100_fcc^2;
end

% FCC�Ŀ�Ѩ��
if weight.vac_energy(1) ~= 0
    evac_fcc_cal = fit_fcc_vac(F, rho, phi, N_fcc_vac, N_fcc, -1*Esub_fcc, r_fcc);    
    out = out + weight.vac_energy(1) * (evac_fcc - evac_fcc_cal)^2/evac_fcc^2;
end

% FCC��100��������죬Bain·��
if weight.fcc_strain(1) ~= 0
    fcc_strain_100_cal = fit_fcc_strain_100(F, rho, phi, fcc_strain_100, a_fcc, Esub_fcc, para.range(4));
    weight_single = weight.fcc_strain(1)/size(fcc_strain_100_cal, 1);
    
    for i = 1:size(fcc_strain_100_cal, 1)
        if fcc_strain_100_cal(i, 2) < fcc_strain_100(i, 2)    % ���������û�����������ڣ������ͷ�
            out = out + weight_single * (fcc_strain_100(i, 2) - fcc_strain_100_cal(i, 2))^2/fcc_strain_100(i, 2)^2;
        elseif fcc_strain_100_cal(i, 2) > fcc_strain_100(i, 3)
            out = out + weight_single * (fcc_strain_100(i, 3) - fcc_strain_100_cal(i, 2))^2/fcc_strain_100(i, 3)^2;
        end
    end
end

% FCC��110���������
if weight.fcc_strain(2) ~= 0
    fcc_strain_110_cal = fit_fcc_strain_110(F, rho, phi, fcc_strain_110, a_fcc, Esub_fcc, para.range(4));
    weight_single = weight.fcc_strain(2)/size(fcc_strain_110_cal, 1);
    
    for i = 1:size(fcc_strain_110_cal, 1)
        if fcc_strain_110_cal(i, 2) < fcc_strain_110(i, 2)    % ���������û�����������ڣ������ͷ�
            out = out + weight_single * (fcc_strain_110(i, 2) - fcc_strain_110_cal(i, 2))^2/fcc_strain_110(i, 2)^2;
        elseif fcc_strain_110_cal(i, 2) > fcc_strain_110(i, 3)
            out = out + weight_single * (fcc_strain_110(i, 3) - fcc_strain_110_cal(i, 2))^2/fcc_strain_110(i, 3)^2;
        end
    end
end

% FCC��111���������
if weight.fcc_strain(3) ~= 0
    fcc_strain_111_cal = fit_fcc_strain_111(F, rho, phi, fcc_strain_111, a_fcc, Esub_fcc, para.range(4));
    weight_single = weight.fcc_strain(3)/size(fcc_strain_111_cal, 1);
    
    for i = 1:size(fcc_strain_111_cal, 1)
        if fcc_strain_111_cal(i, 2) < fcc_strain_111(i, 2)    % ���������û�����������ڣ������ͷ�
            out = out + weight_single * (fcc_strain_111(i, 2) - fcc_strain_111_cal(i, 2))^2/fcc_strain_111(i, 2)^2;
        elseif fcc_strain_111_cal(i, 2) > fcc_strain_111(i, 3)
            out = out + weight_single * (fcc_strain_111(i, 3) - fcc_strain_111_cal(i, 2))^2/fcc_strain_111(i, 3)^2;
        end
    end
end

%% �����ض����͵����ƣ�����---������
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

%% HCP��FCC�Ĺ�ϵ
if weight.prevent_hcp_energy ~= 0
    hcp_energy = fit_prevent_hcp(F, rho, phi, prevent_hcp_energy(1:2), a_fcc, Esub_fcc, para.range(4));
    if hcp_energy < prevent_hcp_energy(3)
        out = out + weight.prevent_hcp_energy * (hcp_energy - prevent_hcp_energy(3))^2/0.01^2;  % ������ѡһ���ϴ���������ж�
    elseif hcp_energy > prevent_hcp_energy(4)
        out = out + weight.prevent_hcp_energy * (hcp_energy - prevent_hcp_energy(4))^2/prevent_hcp_energy(4)^2;
    end
end

%% phi��������״
if weight.phi_well ~= 0
    [diff_cal, max_phi_cal] = fit_phi_well(phi, phi_well);
    if diff_cal < phi_well(4)       % ���ֵ�����˵�һ��������Եڶ������µ����
        out = out + weight.phi_well/2 * abs(diff_cal - phi_well(4));
    end
    if max_phi_cal > phi_well(6)
        out = out + weight.phi_well/2 * abs(max_phi_cal - phi_well(6));
    end
end

%% FCC��Q���Ƶ��
if weight.phonon(1) ~= 0
    frequency_fcc_cal = fit_fcc_frequency(para.DynData, ddphi_fcc, dphi_fcc, ddF_fcc, dF_fcc, ddrho_fcc, drho_fcc, alpha_fcc, r_fcc);
    weight_single = weight.phonon(1)/(size(cons_q, 1)-1);  % ���������������Ȩ��
    out = out + sum(weight_single * (frequency_fcc_cal - cons_q).^2./cons_q.^2);
end
end

function stop = Output(x, optimValues, state)
stop = false;
load Output_lowest para
fval = optimValues.fval;
switch state
    case 'iter'
        switch mod(optimValues.iteration, 2000)      % ÿ1000����һ��IO�������������Լ��ټ����ٶ�
            case 0
                save Output_iter_1 x para fval;
            case 1000
                save Output_iter_2 x para fval;
        end
end
end
