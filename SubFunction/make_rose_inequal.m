%% 根据等式方程，给出不等式方程

function [A, b] = make_rose_inequal(Aeq, beq, A, b, EPS_b, flag_boundary)
A_temp = -1*Aeq;
b_temp = -1*beq + EPS_b;
A = [A; A_temp];
b = [b b_temp];

if flag_boundary == 2
    A_temp = Aeq;
    b_temp = beq + EPS_b;
    A = [A; A_temp];
    b = [b b_temp];
end
end