%% 将约束条件缩并且规整

function F_phi_cons = F_phi_cons_regulate(F_phi_cons)
tol_F = 1e-3;     % F自变量的容忍度
tol_phi = 1e-3;   % phi自变量的容忍度

inequal_coeff = regulate_coeff(F_phi_cons.inequal_coeff, tol_F, tol_phi);       % 先整理不等式
F_phi_cons.inequal_coeff = inequal_coeff;

equal_coeff = regulate_coeff(F_phi_cons.equal_coeff, tol_F, tol_phi);           % 再整理等式
F_phi_cons.equal_coeff = equal_coeff;
end

%% 规整参数
function coeff = regulate_coeff(coeff, tol_F, tol_phi)
MAX = 10000;
idx = find_first_index(coeff);          % 将未使用的数据删除
coeff(idx:end) = [];
record_F = zeros(MAX, 2);
idx_record_F = 1;
record_phi = zeros(MAX, 2);
idx_record_phi = 1;

for i = 1:size(coeff, 1)
    F_coeff = sort_combine(coeff(i).F_coeff, tol_F);
    coeff(i).F_coeff = F_coeff;
    record_F(idx_record_F:idx_record_F+size(F_coeff,1)-1, :) = F_coeff;
    idx_record_F = idx_record_F + size(F_coeff,1);
    
    phi_coeff = sort_combine(coeff(i).phi_coeff, tol_phi);
    coeff(i).phi_coeff = phi_coeff;
    record_phi(idx_record_phi:idx_record_phi+size(phi_coeff,1)-1, :) = phi_coeff;
    idx_record_phi = idx_record_phi + size(phi_coeff,1);
end
record_F(idx_record_F:end, :) = [];
record_phi(idx_record_phi:end, :) = [];

% figure
% hist(record_F(:, 2), 100)
% figure
% hist(record_phi(:, 2), 100)
end

function coeff_new = sort_combine(coeff, tol)
EPS = 1e-6;

[~, idx] = sort(coeff(:, 2));     % 对约束进行排序
coeff = coeff(idx, :);
coeff_new = zeros(size(coeff));

coeff_new(1, :) = coeff(1, :);
idx_new = 1;
idx_old = 1;
if size(coeff, 1) >= 2
    for i = 2:size(coeff, 1)
        if coeff(i, 2) - coeff(i-1, 2) > tol
            coeff_new(idx_new, 1) = sum(coeff(idx_old:i-1, 1));
            coeff_new(idx_new, 2) = mean(coeff(idx_old:i-1, 2));
            idx_new = idx_new + 1;
            idx_old = i;
        end
        if i == size(coeff, 1)
            coeff_new(idx_new, 1) = sum(coeff(idx_old:i, 1));
            coeff_new(idx_new, 2) = mean(coeff(idx_old:i, 2));
        end
    end
end
coeff_new(idx_new+1:end, :) = [];
coeff_new(abs(coeff_new(:, 1)) < EPS, :) = [];       % 如果第一项绝对值过小，也将其删除
end

function idx = find_first_index(coeff)
for i = 1:size(coeff, 1)
    if coeff(i).flag == -1
        idx = i;
        break;
    end
end
end
