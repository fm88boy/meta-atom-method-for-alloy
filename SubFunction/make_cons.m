%% ����Լ��������������Լ������

function [x0, lb, ub, A, b, F_phi_cons] = make_cons(F_data, phi_data, para, cons_data, control, F_phi_cons)
knots = para.knots;     % ��¼�˶���ʽ�Ľڵ����
order = para.order;     % ��¼�˶���ʽ�Ľ״Σ��������������
rho = para.rho;         % ��¼�˵���ܶȺ���
range = para.range;     % ��¼��F��phi������ȡֵ��Χ

alpha_fcc_gsf = para.alpha_fcc_gsf;     % ������FCC�������ƽ��ȱ���ܵļ�������
alpha_fcc = para.alpha_fcc;             % ������FCCƽ��̬�ļ�������
N_fcc = para.N_fcc;
N_fcc_111 = para.N_fcc_111;
N_fcc_110 = para.N_fcc_110;
N_fcc_100 = para.N_fcc_100;

% �Ժ���Ҫ���ⲿ�ִ���ת����ṹ����
drho = (rho(3:end, 2) - rho(1:end-2, 2)) ./ (rho(3:end, 1) - rho(1:end-2, 1));
drho_x = rho(2:end-1, 1);

% ���ڱ���Լ����
a_fcc = cons_data.lattice(1);           % FCC,HCP��BCC�ľ�����
a_hcp = cons_data.lattice(2);
c_hcp = cons_data.lattice(3);
a_bcc = cons_data.lattice(4);
fcc_elastic = cons_data.fcc_elastic;    % FCC,HCP��BCC�ĵ��Գ���
% hcp_elastic = cons_data.hcp_elastic;
% bcc_elastic = cons_data.bcc_elastic;
sub_energy = cons_data.sub_energy;      % FCC,HCP��BCC��������
sub_fcc = sub_energy(1);
sub_hcp = sub_energy(2);
sub_bcc = sub_energy(3);

% ���ڱ���Լ������
F_phi_range = control.F_phi_range;
fcc_relation_C12_C44 = control.fcc_relation_C12_C44;
dF_restrict = control.dF_restrict;
fcc_dF_restrict_equal = control.fcc_dF_restrict_equal;
phi_restrict = control.phi_restrict;
dphi_restrict = control.dphi_restrict;
ddphi_restrict = control.ddphi_restrict;
fcc_rose_equation = control.fcc_rose_equation;
fcc_100_tension = control.fcc_100_tension;
hcp_rose_equation = control.hcp_rose_equation;
bcc_rose_equation = control.bcc_rose_equation;
fcc_sf_relation = control.fcc_sf_relation;
fcc_surface_relation = control.fcc_surface_relation;

x0 = [F_data, phi_data];        % �Ż���ʼ��

% ����FCC��һЩ��������
radius_fcc = alpha_fcc*a_fcc;               % FCC�İ뾶
Aeq_phi_fcc = zeros(size(phi_data));        % FCCƽ����phi����
for i = 1:size(N_fcc, 2)
    Aeq_phi_fcc = Aeq_phi_fcc + 0.5*N_fcc(i)*cal_phi(radius_fcc(i), knots, order);
end
rho_fcc = sum(N_fcc .* interp1(rho(:,1), rho(:,2), radius_fcc));   % FCCƽ���ĵ���ܶ�
drho_fcc = interp1(drho_x, drho, radius_fcc);

% ����F������phi���������޺�����
lb_F = -1*F_phi_range(1)*ones(size(F_data));
ub_F = F_phi_range(1)*ones(size(F_data));
lb_phi = -1*F_phi_range(2)*ones(size(phi_data));
ub_phi = F_phi_range(2)*ones(size(phi_data));
lb = [lb_F, lb_phi];     % lb�Ļ�������F��-1������phi��-100
ub = [ub_F, ub_phi];     % ub�Ļ�������F��1������phi��100

A = [];     % ��ʼ��
b = [];

% ����F�����Ĳ���ʽԼ��,�˴�Ҫ��ʵ���C12-C44<0.5(C12_cons-C44_cons)�Ĳ��
% �����ԣ���һ������������C12<C44�����
if fcc_relation_C12_C44(1) == 1
    xi = rho_fcc;
    delta = [2, 2, 12, 8, 20];
    omega = a_fcc^3/4;
    V = sum(drho_fcc .* delta ./ alpha_fcc * a_fcc);
    A_temp(1, :) = [cal_ddF(xi), zeros(size(phi_data))];
    b_temp = fcc_relation_C12_C44(2) * (fcc_elastic(2)-fcc_elastic(3)) * omega/V^2;
    A = [A; A_temp];
    b = [b, b_temp];
end

% ��֤F������[dF_restrict(3),dF_restrict(4)]����������һ�׵���С��0
if dF_restrict(1) == 1
    N_A = dF_restrict(2);
    xi = linspace(dF_restrict(3), dF_restrict(4), N_A);
    A_temp = zeros(N_A, size(x0, 2));
    for i = 1:N_A
        A_temp(i, :) = [cal_dF(xi(i)), zeros(size(phi_data))];
    end
    b_temp = 0 * ones(size(xi));        % ��Ϊ��С��0��,�������Ϊ0*ones
    A = [A; A_temp];
    b = [b, b_temp];
end

% �Ƿ�����dF_rho_fccֵ�����Ƶ����ֵ�������Ĳ�������
if fcc_dF_restrict_equal(1) == 1
    dF_cons = fcc_dF_restrict_equal(2);
    N_A = fcc_dF_restrict_equal(3);
    xi = linspace(rho_fcc-1, rho_fcc+1, N_A);
    A_temp = zeros(N_A, size(x0, 2));
    for i = 1:N_A
        A_temp(i, :) = [cal_dF(xi(i)), zeros(size(phi_data))];
    end
    b_temp = dF_cons * ones(size(xi));
    A = [A; A_temp];
    b = [b, b_temp];
end

% ����phi�����ڵ�һ��knots��������ֵ
if phi_restrict(1) == 1
    A_temp = [zeros(size(F_data)), cal_phi(knots(1), knots, order)];
    b_temp = phi_restrict(2);
    F_phi_cons = F_phi_cons_record([0, 0], [1, knots(1)], [1, phi_restrict(2)], F_phi_cons);
    A = [A; A_temp];
    b = [b b_temp];
end


% ����phi�����ڲ�ͬλ�ô���һ�׵���С��0
if dphi_restrict(1) == 1
    restrict_point = dphi_restrict(2:end);
    A_temp = zeros(size(restrict_point, 2), size(x0, 2));
    b_temp = zeros(1, size(restrict_point, 2));
    for i = 1:size(restrict_point, 2)
        A_temp(i, :) = [zeros(size(F_data)), cal_dphi(restrict_point(i), knots, order)];
        b_temp(i) = 0;
    end
    A = [A; A_temp];
    b = [b b_temp];
end

% ��һ���ķ�Χ�ڣ����׵������ܳ���һ���޶ȣ�ddphi_restrict(2)
if ddphi_restrict(1) == 1
    N_A = ddphi_restrict(3);
    xi = linspace(ddphi_restrict(4), ddphi_restrict(5), N_A);
    A_temp = zeros(N_A, size(x0, 2));
    b_temp = ddphi_restrict(2) * ones(1, N_A);
    for i = 1:N_A
        A_temp(i, :) = [zeros(size(F_data)), cal_ddphi(xi(i), knots, order)];
    end
    A = [A; A_temp];
    b = [b b_temp];
    
    A_temp = zeros(N_A, size(x0, 2));
    b_temp = ddphi_restrict(2) * ones(1, N_A);
    for i = 1:N_A
        A_temp(i, :) = [zeros(size(F_data)), -1*cal_ddphi(xi(i), knots, order)];
    end
    A = [A; A_temp];
    b = [b b_temp];
end

% ������ʽԼ����������FCC����������ϣ�rose����
if fcc_rose_equation(1) == 1
    N_equ = fcc_rose_equation(2);     % ѡ�񼸸�������Լ������
    scan_span = [-fcc_rose_equation(3), fcc_rose_equation(3)];      % ɨ�������
    r_wse = power(3/16/pi, 1/3) * a_fcc;         % ������ȼ����
    K = (fcc_elastic(1) + 2*fcc_elastic(2))/3;   % ����C11��C12��������ģ��
    para_l = sqrt(sub_fcc/(12*pi*K*r_wse));     % ����������(��ֵ)�;�������������м����l
    r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);
    a_norm = linspace(scan_span(1), scan_span(2), N_equ);
    energy_ref = sub_fcc * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);
    
    Aeq = zeros(N_equ, size(x0, 2));
    beq = zeros(1, N_equ);
    
    xi = linspace(-0.7, 0.7, N_equ);
    EPS_b = fcc_rose_equation(4)*(exp(6*abs(xi))/exp(6)*500+1);     % ��������ʽ��������
    for i = 1 : N_equ
        a_temp = power(16*pi/3, 1/3) * r(i);    % ������FCC����
        radius_this = alpha_fcc*a_temp;
        rho_this = sum(N_fcc.*interp1(rho(:,1), rho(:,2), radius_this));   % �������Ӧ�ĵ���ܶ�
        beq(i) = energy_ref(i);
        Aeq_phi = zeros(size(phi_data));
        for j = 1:size(radius_this, 2)
            Aeq_phi = Aeq_phi + 0.5*N_fcc(j)*cal_phi(radius_this(j), knots, order);
        end
        Aeq(i, :) = [cal_F(rho_this), Aeq_phi];

        F_coeff = [1, rho_this];
        phi_coeff = [0.5*N_fcc', radius_this'];
        b_coeff = [1, beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % ����Լ��
        F_coeff = [-1, rho_this];
        phi_coeff = [-0.5*N_fcc', radius_this'];
        b_coeff = [1, -beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    end
    [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, fcc_rose_equation(5));   % ����һ������Rose���ߵĲ���ʽ    
end

% ������ʽԼ����������HCP����������ϣ�rose���̡�ȫ�̱���a�仯����c/a����
if hcp_rose_equation(1) == 1
    N_equ = hcp_rose_equation(2);     % ѡ�񼸸�������Լ������
    scan_span = [-hcp_rose_equation(3), hcp_rose_equation(3)];      % ɨ�������
    r_wse = power(3*sqrt(3)/(32*pi)*(c_hcp/a_hcp), 1/3) * a_hcp;         % ������ȼ����
    % K = (fcc_elastic(1) + 2*fcc_elastic(2))/3;   % ����C11��C12��������ģ��
    % ��������Ҫʹ��FCC������ģ������Ϊ��FCCΪĸ�࣬Ϊ�����鷳��ֱ�Ӳ���FCC�����ģ��
    
    para_l = sqrt(sub_hcp/(12*pi*K*r_wse));     % ����������(��ֵ)�;�������������м����l
    r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);      % WS�뾶��ɨ�跶Χ
    a_norm = linspace(scan_span(1), scan_span(2), N_equ);
    energy_ref = sub_hcp * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);    % Ŀ�������Ĳο�ֵ
    
    Aeq = zeros(N_equ, size(x0, 2));
    beq = zeros(1, N_equ);
    
    xi = linspace(-0.7, 0.7, N_equ);
    EPS_b = hcp_rose_equation(4)*(exp(6*abs(xi))/exp(6)*500+1);
    for i = 1 : N_equ
        a_temp = 1/power(3*sqrt(3)/(32*pi)*(c_hcp/a_hcp), 1/3) * r(i);    % ������HCP����
        c_temp = c_hcp/a_hcp * a_temp;
        [radius_this, N_this] = cal_hcp_radius(a_temp, c_temp, range(4));         % ����µ�HCP�뾶��ԭ����
        rho_this = sum(N_this.*interp1(rho(:,1), rho(:,2), radius_this));         % �������Ӧ�ĵ���ܶ�
        beq(i) = energy_ref(i);
        Aeq_phi = zeros(size(phi_data));
        for j = 1:size(radius_this, 2)
            Aeq_phi = Aeq_phi + 0.5*N_this(j)*cal_phi(radius_this(j), knots, order);
        end
        Aeq(i, :) = [cal_F(rho_this), Aeq_phi];
        
        F_coeff = [1, rho_this];
        phi_coeff = [0.5*N_this', radius_this'];
        b_coeff = [1, beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % ����Լ��
        F_coeff = [-1, rho_this];
        phi_coeff = [-0.5*N_this', radius_this'];
        b_coeff = [1, -beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    end
    [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, hcp_rose_equation(5));   % ����һ������Rose���ߵĲ���ʽ
end

% ������ʽԼ����������BCC����������ϣ�rose����
if bcc_rose_equation(1) == 1
    N_equ = bcc_rose_equation(2);     % ѡ�񼸸�������Լ������
    scan_span = [-bcc_rose_equation(3), bcc_rose_equation(3)];      % ɨ�������
    r_wse = power(3/8/pi, 1/3) * a_bcc;         % ������ȼ����
    % K = (fcc_elastic(1) + 2*fcc_elastic(2))/3;   % ����C11��C12��������ģ��
    % ������Ҳ��ʹ����FCC������ģ������������������Fe50Mn30Cr10Co10��HEA�Ͻ�
    
    para_l = sqrt(sub_bcc/(12*pi*K*r_wse));     % ����������(��ֵ)�;�������������м����l
    r = linspace(r_wse+scan_span(1)*para_l, r_wse+scan_span(2)*para_l, N_equ);
    a_norm = linspace(scan_span(1), scan_span(2), N_equ);
    energy_ref = sub_bcc * (-1-a_norm-0.05*a_norm.^3) .* exp(-a_norm);
    
    Aeq = zeros(N_equ, size(x0, 2));
    beq = zeros(1, N_equ);
    
    xi = linspace(-0.7, 0.7, N_equ);
    EPS_b = bcc_rose_equation(4)*(exp(6*abs(xi))/exp(6)*500+1);
    for i = 1 : N_equ
        a_temp = power(8*pi/3, 1/3) * r(i);    % ������BCC����
        [radius_this, N_this] = cal_bcc_radius(a_temp, range(4));      % �������Ӧ�ļ�������
        rho_this = sum(N_this .* interp1(rho(:,1), rho(:,2), radius_this));   % �������Ӧ�ĵ���ܶ�
        beq(i) = energy_ref(i);
        Aeq_phi = zeros(size(phi_data));
        for j = 1:size(radius_this, 2)
            Aeq_phi = Aeq_phi + 0.5*N_this(j)*cal_phi(radius_this(j), knots, order);
        end
        Aeq(i, :) = [cal_F(rho_this), Aeq_phi];

        F_coeff = [1, rho_this];
        phi_coeff = [0.5*N_this', radius_this'];
        b_coeff = [1, beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);    % ����Լ��
        F_coeff = [-1, rho_this];
        phi_coeff = [-0.5*N_this', radius_this'];
        b_coeff = [1, -beq(i)+EPS_b(i)];
        F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons);
    end
    [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, bcc_rose_equation(5));
end

% ���������FCC����100���������Լ��
if fcc_100_tension(1) == 1
    N_inequ = (size(fcc_100_tension, 2)-1)/2;
    index = 2;      % �����ָ��
    for i = 1 : N_inequ
        [A, b, F_phi_cons] = make_fcc_strain_100_inequal(A, b, rho, phi_data, Aeq_phi_fcc, rho_fcc, knots, order, a_fcc, fcc_100_tension(index:index+1), ...
            fcc_elastic(1), range(4), N_fcc, radius_fcc, F_phi_cons);
        index = index + 2;
    end
end

% ����FCC����Ĳ���ʽԼ��������Ҫ�����ƽ��ȱ��������һ���������ʽ
% Egsf_0_1 < Egsf_0_2 < Egsf_0_3 < Egsf_0_4 < Egsf_0_5(���ȶ������)
% Egsf_1_0(�����) < Egsf_0_5(���ȶ������)

if fcc_sf_relation(1) == 1
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [1, 2], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [2, 3], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [3, 4], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [4, 5], F_phi_cons);
    [A, b, F_phi_cons] = make_fcc_sf_inequal(A, b, rho, phi_data, knots, order, a_fcc, alpha_fcc_gsf, [6, 5], F_phi_cons);
end

% ��������FCC�����������Լ����Ҫ��E111<E100<E110
if fcc_surface_relation(1) == 1
    [A, b, F_phi_cons] = make_fcc_surface_inequal(A, b, rho, phi_data, Aeq_phi_fcc, rho_fcc, knots, order, a_fcc, N_fcc_111, ...
        N_fcc_110, N_fcc_100, N_fcc, radius_fcc, F_phi_cons);
end

end
