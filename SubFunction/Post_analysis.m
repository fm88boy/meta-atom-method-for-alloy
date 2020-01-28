%% �������

function [out_txt, para_out] = Post_analysis(data, para, is_write, write_para)
close all

%% ���´����ֱ�ӴӶ�Ӧ��Inner_fun����
% D171027
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
    out = out + weight.sf_energy(1) * 0.6 * (esf_fcc - esf_fcc_cal)^2/esf_fcc^2;        % ���ڲ���ܣ�����60%��Ȩ��
    out = out + weight.sf_energy(1) * 0.4 * (eusf_fcc - eusf_fcc_cal)^2/eusf_fcc^2;     % ���ڲ��ȶ�����ܣ�����40%��Ȩ��
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

%% ��ʼ����
EPS = 1e-4;
if norm(rho(:, 1) - phi(:, 1)) < EPS        % ���rho������phi��������ƥ�䣬��Ҫ��������rho����
    rho_temp = rho;
    rho = zeros(size(phi));
    rho(:, 1) = phi(:, 1);
    rho(:, 2) = interp1(rho_temp(:, 1), rho_temp(:, 2), rho(:, 1));
    nan_idx_1 = isnan(rho(:, 2)) & rho(:, 1) < 2;       % ��ǰ���nan�����
    rho(nan_idx_1, 2) = max(rho_temp(:, 2));
    nan_idx_2 = isnan(rho(:, 2)) & rho(:, 1) > 4;
    rho(nan_idx_2, 2) = 0;
end

if is_write == 1
    % ��ʼ���ƻ���������F������phi����
    figure_idx = figure('InvertHardcopy','off','Color',[1 1 1]);
    axes_idx = axes('Parent', figure_idx, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 14);
    box(axes_idx, 'on');
    hold(axes_idx, 'on');
    plot(F(:, 1), F(:, 2), 'r', 'LineWidth', 2)
    xlim([para.range(1), 40])
    ylabel('F', 'FontWeight', 'bold', 'FontSize', 15.4)
    xlabel('rho', 'FontWeight', 'bold', 'FontSize', 15.4)
    
    figure_idx = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
    axes_idx = axes('Parent', figure_idx, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 14);
    box(axes_idx, 'on');
    hold(axes_idx, 'on');
    plot(phi(:, 1), phi(:, 2), 'r', 'LineWidth', 2)
    xlim([2, para.range(4)])
    ylim([-1 2])
    ylabel('phi','FontWeight','bold','FontSize',15.4)
    xlabel('radius','FontWeight','bold','FontSize',15.4)
end

%% fcc��100��������
if weight.fcc_strain(1) ~= 0 && is_write == 1
    plot_fcc_strain_100 = zeros(40, 2);          % ��������������ߵ�
    plot_fcc_strain_100(:, 1) = linspace(-0.1, 0.15, 40)';
    plot_fcc_strain_100_data = fit_fcc_strain_100(F, rho, phi, plot_fcc_strain_100, a_fcc, Esub_fcc, para.range(4));
    
    figure_idx = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
    axes_idx = axes('Parent', figure_idx, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 14);
    box(axes_idx, 'on');
    hold(axes_idx, 'on');
    plot(plot_fcc_strain_100_data(:, 1), plot_fcc_strain_100_data(:, 2), 'LineWidth', 2)
    xlabel('e_{100}', 'FontWeight', 'bold', 'FontSize', 15.4)
    ylabel('energy', 'FontWeight', 'bold', 'FontSize', 15.4)
end

%% fcc��110��������
if weight.fcc_strain(2) ~= 0 && is_write == 1
    plot_fcc_strain_110 = zeros(40, 2);          % ��������������ߵ�
    plot_fcc_strain_110(:, 1) = linspace(-0.1, 0.15, 40)';
    plot_fcc_strain_110_data = fit_fcc_strain_110(F, rho, phi, plot_fcc_strain_110, a_fcc, Esub_fcc, para.range(4));
    
    figure_idx = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
    axes_idx = axes('Parent', figure_idx, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 14);
    box(axes_idx, 'on');
    hold(axes_idx, 'on');
    plot(plot_fcc_strain_110_data(:, 1), plot_fcc_strain_110_data(:, 2), 'LineWidth', 2)
    xlabel('e_{110}', 'FontWeight', 'bold', 'FontSize', 15.4)
    ylabel('energy', 'FontWeight', 'bold', 'FontSize', 15.4)
end

%% fcc��111��������
if weight.fcc_strain(3) ~= 0 && is_write == 1
    plot_fcc_strain_111 = zeros(40, 2);          % ��������������ߵ�
    plot_fcc_strain_111(:, 1) = linspace(-0.1, 0.15, 40)';
    plot_fcc_strain_111_data = fit_fcc_strain_111(F, rho, phi, plot_fcc_strain_111, a_fcc, Esub_fcc, para.range(4));
    
    figure_idx = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
    axes_idx = axes('Parent', figure_idx, 'LineWidth', 2, 'FontWeight', 'bold', 'FontSize', 14);
    box(axes_idx, 'on');
    hold(axes_idx, 'on');
    plot(plot_fcc_strain_111_data(:, 1), plot_fcc_strain_111_data(:, 2), 'LineWidth', 2)
    xlabel('e_{111}', 'FontWeight', 'bold', 'FontSize', 15.4)
    ylabel('energy', 'FontWeight', 'bold', 'FontSize', 15.4)
end

%% FCC��Rose����
if a_fcc ~= 0 && is_write == 1
    [a_real, energy_ref, energy] = fit_fcc_rose(F, rho, phi, 11, [-0.15, 0.15], a_fcc, Esub_fcc, C11_fcc, C12_fcc, alpha_fcc, N_fcc);
    
    figure_idx = figure('InvertHardcopy','off','Color',[1 1 1]);        % ���Ƴ�Rose����
    axes_idx = axes('Parent',figure_idx,'LineWidth',2,'FontWeight','bold', 'FontSize',14);
    box(axes_idx,'on');
    hold(axes_idx,'on');
    plot(a_real, energy_ref, 'b:', a_real, energy, 'r', a_real, energy, 'o', 'LineWidth',2)
    ylabel('Energy', 'FontWeight','bold','FontSize',15.4)
    xlabel('Lattice constant', 'FontWeight','bold','FontSize',15.4)
end

%% HCP��Rose����
if a_hcp ~= 0 && is_write == 1
    K_hcp = (C11_fcc + 2*C12_fcc)/3;    % ���棬�ⲿ�ִ���ֻ��������Fe50Mn30Cr10Co10��ʹ��
    [a_real, energy_ref, energy] = fit_hcp_rose(F, rho, phi, 11, [-0.15, 0.15], a_hcp, c_hcp, Esub_hcp, K_hcp, para.range(4));
    
    figure_idx = figure('InvertHardcopy','off','Color',[1 1 1]);        % ���Ƴ�Rose����
    axes_idx = axes('Parent',figure_idx,'LineWidth',2,'FontWeight','bold', 'FontSize',14);
    box(axes_idx,'on');
    hold(axes_idx,'on');
    plot(a_real, energy_ref, 'b:', a_real, energy, 'r', a_real, energy, 'o', 'LineWidth',2)
    ylabel('Energy', 'FontWeight','bold','FontSize',15.4)
    xlabel('Lattice constant', 'FontWeight','bold','FontSize',15.4)
end

%% BCC��Rose����
if a_bcc ~= 0 && is_write == 1
    K_bcc = (C11_fcc + 2*C12_fcc)/3;    % ���棬�ⲿ�ִ���ֻ��������Fe50Mn30Cr10Co10��ʹ��
    [a_real, energy_ref, energy] = fit_bcc_rose(F, rho, phi, 11, [-0.15, 0.15], a_bcc, Esub_bcc, K_bcc, para.range(4));
    
    figure_idx = figure('InvertHardcopy','off','Color',[1 1 1]);        % ���Ƴ�Rose����
    axes_idx = axes('Parent',figure_idx,'LineWidth',2,'FontWeight','bold', 'FontSize',14);
    box(axes_idx,'on');
    hold(axes_idx,'on');
    plot(a_real, energy_ref, 'b:', a_real, energy, 'r', a_real, energy, 'o', 'LineWidth',2)
    ylabel('Energy', 'FontWeight','bold','FontSize',15.4)
    xlabel('Lattice constant', 'FontWeight','bold','FontSize',15.4)
end

%% ������ⲿ�ִ������������

% ���ȰѼ������ݵĽṹ�廯
cal_data = para.cons_data;
cal_data.lattice = [a_fcc_cal, [a_hcp_cal, c_hcp_cal], a_bcc_cal];          % ������������ΪFCC��HCP_a��HCP_c��BCC
cal_data.fcc_elastic = [C11_fcc_cal, C12_fcc_cal, C44_fcc_cal];             % FCC���Գ���������˳������ΪC11,C12,C44
cal_data.sub_energy = abs([fcc_energy_cal, hcp_energy_cal, bcc_energy_cal]);% ������
cal_data.sf_energy = [esf_fcc_cal, 0, 0];                                   % �����
cal_data.usf_energy = [eusf_fcc_cal, 0, 0];                                 % ���ȶ������
cal_data.fcc_surface_energy = [e111_fcc_cal, e110_fcc_cal, e100_fcc_cal];   % FCC�����ܣ�����Ϊ111,110��100��
cal_data.vac_energy = [evac_fcc_cal, 0, 0];                                 % ��Ѩ�γ���
cal_data.fcc_strain_100 = fcc_strain_100_cal;                               % FCC��100·����Ӧ�䣬Ŀ��������eV��
cal_data.fcc_strain_110 = fcc_strain_110_cal;                               % FCC��110·����Ӧ�䣬Ŀ��������eV��
cal_data.fcc_strain_111 = fcc_strain_111_cal;                               % FCC��111·����Ӧ�䣬Ŀ��������eV��
cal_data.fcc_specific_energy = fcc_specific_energy_cal;                     % FCC�����⹹�͵�����
cal_data.phi_well = [diff_cal, max_phi_cal];                                % phi������״�����ã�����1������2����С���,������phiֵ

para_out.out = out;
para_out.cons_data = para.cons_data;
para_out.cons_data.fcc_specific_energy = para.fcc_specific_energy;          % �����⹹�͵������ٻ�����
para_out.cons_data.fcc_elastic = para_out.cons_data.fcc_elastic*160;
para_out.cons_data.sf_energy = para_out.cons_data.sf_energy*16*1e3;
para_out.cons_data.usf_energy = para_out.cons_data.usf_energy*16*1e3;
para_out.cons_data.fcc_surface_energy = para_out.cons_data.fcc_surface_energy*16*1e3;

para_out.cal_data = cal_data;
para_out.cal_data.fcc_elastic = para_out.cal_data.fcc_elastic*160;
para_out.cal_data.sf_energy = para_out.cal_data.sf_energy*16*1e3;
para_out.cal_data.usf_energy = para_out.cal_data.usf_energy*16*1e3;
para_out.cal_data.fcc_surface_energy = para_out.cal_data.fcc_surface_energy*16*1e3;

text_a = sprintf('Test results:\na_fcc=%.5f  C11_fcc=%.2f  C12_fcc=%.2f  C44_fcc=%.2f  ', ...
    a_fcc_cal, [C11_fcc_cal, C12_fcc_cal, C44_fcc_cal]*160);
text_b = sprintf('esub_fcc=%.5f  esf_fcc=%.2f  eusf_fcc=%.2f  e111_fcc=%.2f  e110_fcc=%.2f  e100_fcc=%.2f  evac_fcc=%.2f\n', ...
    abs(fcc_energy_cal), [esf_fcc_cal, eusf_fcc_cal, e111_fcc_cal, e110_fcc_cal, e100_fcc_cal]*16*1e3, evac_fcc_cal);
text_c = sprintf('a_hcp=%.5f  c_hcp=%.5f  esub_hcp=%.5f  a_bcc=%.5f  esub_bcc=%.5f\n', ...
    a_hcp_cal, c_hcp_cal, abs(hcp_energy_cal), a_bcc_cal, abs(bcc_energy_cal));

text_d = sprintf('Cons results:\na_fcc=%.5f  C11_fcc=%.2f  C12_fcc=%.2f  C44_fcc=%.2f  ', ...
    a_fcc, [C11_fcc, C12_fcc, C44_fcc]*160);
text_e = sprintf('esub_fcc=%.5f  esf_fcc=%.2f  eusf_fcc=%.2f  e111_fcc=%.2f  e110_fcc=%.2f  e100_fcc=%.2f  evac_fcc=%.2f\n', ...
    Esub_fcc, [esf_fcc, eusf_fcc, e111_fcc, e110_fcc, e100_fcc]*16*1e3, evac_fcc);
text_f = sprintf('a_hcp=%.5f  c_hcp=%.5f  esub_hcp=%.5f  a_bcc=%.5f  esub_bcc=%.5f\n', ...
    a_hcp, c_hcp, Esub_hcp, a_bcc, Esub_bcc);

temp_1 = sprintf('e_fcc_100=%.2f\nTest: energy=%.2f   Cons: energy_low=%.2f   energy_high=%.2f\n', ...
    [fcc_strain_100_cal, fcc_strain_100(:, 2:3)]');
temp_2 = sprintf('e_fcc_110=%.2f\nTest: energy=%.2f   Cons: energy_low=%.2f   energy_high=%.2f\n', ...
    [fcc_strain_110_cal, fcc_strain_110(:, 2:3)]');
temp_3 = sprintf('e_fcc_111=%.2f\nTest: energy=%.2f   Cons: energy_low=%.2f   energy_high=%.2f\n', ...
    [fcc_strain_111_cal, fcc_strain_111(:, 2:3)]');
text_g = sprintf('\nfcc_100 path:\n%sfcc_110 path:\n%sfcc_111 path:\n%s\n', temp_1, temp_2, temp_3);

index_id = (1:1:size(fcc_specific_energy_cal, 1))';
temp = sprintf('fcc_specific_config_id=%d   Test: %.4f   Cons: energy_low=%.4f   energy_high=%.4f\n', ...
    [index_id, fcc_specific_energy_cal, fcc_specific_energy]');
text_g = sprintf('%s%s\n', text_g, temp);

text_h = sprintf('Well:\nTest: diff_energy= %.2f   max_phi = %.2f\nCons: diff_energy= %.2f   max_phi = %.2f\n', ...
    diff_cal, max_phi_cal, phi_well(4), phi_well(6));

out_txt = sprintf('%s%s%s%s%s%s%s%s\n', text_a, text_b, text_c, text_d, text_e, text_f, text_g, text_h);    % ����ı�
% disp(out_txt)

if is_write == 1
    filename = sprintf('%s\\Work', write_para.path);
    if exist(filename, 'dir') ~= 7
        mkdir(filename);
    end
    figure(1)
    filename = sprintf('%s\\Work\\%d-F.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(2)
    filename = sprintf('%s\\Work\\%d-phi.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(3)
    filename = sprintf('%s\\Work\\%d-fcc-100.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(4)
    filename = sprintf('%s\\Work\\%d-fcc-110.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(5)
    filename = sprintf('%s\\Work\\%d-fcc-111.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(6)
    filename = sprintf('%s\\Work\\%d-rose-fcc.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(7)
    filename = sprintf('%s\\Work\\%d-rose-hcp.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    figure(8)
    filename = sprintf('%s\\Work\\%d-rose-bcc.jpg', write_para.path, write_para.write_index);  % ֻ������Windows
    print('-djpeg', '-opengl', '-r200', filename)
    filename = sprintf('%s\\Work\\%d-inform.txt', write_para.path, write_para.write_index);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', out_txt);
    fclose(fid);
    
    filename = sprintf('%s\\Work\\%d-HEA.eam.fs', write_para.path, write_para.write_index);
    Write_EAM(F, rho, phi, {write_para.element{1}, write_para.element{2}}, [write_para.element{3}, write_para.element{4}, a_fcc], ...
        [size(F, 1), F(2, 1), size(phi, 1), phi(2, 1), phi(end, 1)], filename, write_para.comment)
end
end
