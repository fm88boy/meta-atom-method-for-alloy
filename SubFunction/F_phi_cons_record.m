%% 根据约束条件，记录F和phi的所有约束条件

function F_phi_cons = F_phi_cons_record(F_coeff, phi_coeff, b_coeff, F_phi_cons)
% F_coeff：第一列F值的系数，第二列F值的自变量
% phi_coeff：第一列phi值的系数，第二列phi值的自变量
% b_coeff：第一项方程类型（1：不等式(<)；2：等式），第二项方程的值

if b_coeff(1) == 1
    inequal_coeff = F_phi_cons.inequal_coeff;
    idx = find_first_index(inequal_coeff);          % 找到第一个空数据
    inequal_coeff(idx).F_coeff = F_coeff;
    inequal_coeff(idx).phi_coeff = phi_coeff;
    inequal_coeff(idx).b_coeff = b_coeff(2);
    inequal_coeff(idx).flag = 1;        % 标记为已保存
    F_phi_cons.inequal_coeff = inequal_coeff;
elseif b_coeff(1) == 2
    equal_coeff = F_phi_cons.equal_coeff;
    idx = find_first_index(equal_coeff);          % 找到第一个空数据
    equal_coeff(idx).F_coeff = F_coeff;
    equal_coeff(idx).phi_coeff = phi_coeff;
    equal_coeff(idx).b_coeff = b_coeff(2);
    equal_coeff(idx).flag = 1;        % 设置为已保存
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