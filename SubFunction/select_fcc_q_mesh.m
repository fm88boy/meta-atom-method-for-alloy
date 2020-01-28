%% ����Լ����Ŀ�����������ݿ��ж�ȡ��Լ������

function [D, cons_q] = select_fcc_q_mesh(cons_q_mesh, q_list, D_list)
D = cell(size(cons_q_mesh, 1), 1);
cons_q = zeros(size(cons_q_mesh, 1), 1);

EPS = 1e-4;
for i = 1:size(cons_q_mesh, 1)
    for j = 1:size(q_list, 1)
        if norm(cons_q_mesh(i, 1:3) - q_list(j, :)) < EPS
            D{i} = D_list(:, :, j);     % �����ݿ��а���Ϣ����ȡ����
            break;
        end
    end
    cons_q(i) = cons_q_mesh(i, 4);
end
end
