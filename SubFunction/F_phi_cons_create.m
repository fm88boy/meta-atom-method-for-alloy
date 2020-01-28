%% ���������Լ��������������һ��struct������

function F_phi_cons = F_phi_cons_create(max_num)
% max_num����������������Լ����
F_coeff = cell(max_num, 1);
phi_coeff = cell(max_num, 1);
b_coeff = cell(max_num, 1);
inequal_coeff = struct('flag', -1, 'F_coeff', F_coeff, 'phi_coeff', phi_coeff, 'b_coeff', b_coeff);
equal_coeff = struct('flag', -1, 'F_coeff', F_coeff, 'phi_coeff', phi_coeff, 'b_coeff', b_coeff);

F_phi_cons.inequal_coeff = inequal_coeff;
F_phi_cons.equal_coeff = equal_coeff;
end