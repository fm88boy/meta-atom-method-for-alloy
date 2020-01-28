%% ����111�᷽����м���

function [e_expan, energy] = loading_fcc_strain_111(F, rho, phi, strain, a, r_cut, alpha_fcc_strain_111, N_fcc_strain_111, length)
% ��������һ��������Ϊa��FCC��������ĳ������Ӧ��Ϊstrainʱϵͳ�������������������������ɷ���
bnd = [-1.5*abs(strain), 0.5*abs(strain)];    % �������Ϊ-1.5*e_zz��0.5*e_zz��ʵ���ϣ�����ֵӦ����-0.5*e_zz
options = optimset('Display', 'off');
[e_expan, energy] = fminbnd(@(e_expan)expansion_fcc_111(F, rho, phi, a, e_expan, strain, r_cut, alpha_fcc_strain_111, N_fcc_strain_111, length), ...
    bnd(1), bnd(2), options);
end

function energy = expansion_fcc_111(F, rho, phi, a_fcc, e_expan, e_zz, r_cut, alpha_fcc_strain_111, N_fcc_strain_111, length)
a_fcc_cal = a_fcc * (1+e_expan);
h_cal = sqrt(3)/3*a_fcc * (1+e_zz);     % ���������Ҫ���ݷ�����ı�

radius = zeros(1, length);
N = zeros(1, length);
input_idx = 1;
for i = 1:size(alpha_fcc_strain_111, 2)
    this_length = size(alpha_fcc_strain_111{i}, 2);     % �����
    radius(input_idx:input_idx+this_length-1) = sqrt((alpha_fcc_strain_111{i}*a_fcc_cal).^2 + ((i-1)*h_cal)^2);
    N(input_idx:input_idx+this_length-1) = N_fcc_strain_111{i};
    input_idx = input_idx + this_length;
end

idx = radius < r_cut;
radius = radius(idx);
N = N(idx);

phi_fcc = sum(N .* interp1(phi(:,1), phi(:,2), radius));
rho_fcc = sum(N .* interp1(rho(:,1), rho(:,2), radius));
energy = interp1(F(:, 1), F(:, 2), rho_fcc) + 0.5*phi_fcc;
end
