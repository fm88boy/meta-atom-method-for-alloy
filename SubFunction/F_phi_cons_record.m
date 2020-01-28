%% ����Լ����������¼F��phi������Լ������

function F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons)
% F_coeff����һ��Fֵ��ϵ�����ڶ���Fֵ���Ա���
% phi_coeff����һ��phiֵ��ϵ�����ڶ���phiֵ���Ա���
% b_coeff����һ������ͣ�1������ʽ(<)��2����ʽ�����ڶ���̵�ֵ

if b_coeff(1) == 1
    inequal_coeff = F_phi_cons.inequal_coeff;
    idx = find_first_index(inequal_coeff);          % �ҵ���һ��������
    inequal_coeff(idx).F_coeff = F_coeff;
    inequal_coeff(idx).phi_coeff = phi_coeff;
    inequal_coeff(idx).b_coeff = b_coeff(2);
    inequal_coeff(idx).flag = 1;        % ���Ϊ�ѱ���
    F_phi_cons.inequal_coeff = inequal_coeff;
elseif b_coeff(1) == 2
    equal_coeff = F_phi_cons.equal_coeff;
    idx = find_first_index(equal_coeff);          % �ҵ���һ��������
    equal_coeff(idx).F_coeff = F_coeff;
    equal_coeff(idx).phi_coeff = phi_coeff;
    equal_coeff(idx).b_coeff = b_coeff(2);
    equal_coeff(idx).flag = 1;        % ����Ϊ�ѱ���
    F_phi_cons.equal_coeff = equal_coeff;
end
end

function idx = find_first_index(coeff)
for i = 1:size(coeff, 1)
    if coeff(i).flag == -1
        idx = i;
        break;
    end
end
end