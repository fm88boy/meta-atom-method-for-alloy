% 根据半径、节点和阶次，计算出phi矩阵
function out = cal_phi(radius, knots, order)
% radius：原子半径；knots：节点；order：多项式阶次
num_order = order(2) - order(1) + 1;    % 对于一个节点，自变量的个数
out = zeros(1, (size(knots, 2)-1)*num_order);         % 第一个phi的knots是用来确定下限的

for i = 1:size(out, 2)
    this_knots = ceil(i/num_order)+1;   % 对应的节点
    this_order = order(1) + mod(i-1, num_order);    % 对应的阶次
    out(i) = heaviside(knots(this_knots)-radius)*(knots(this_knots)-radius).^this_order;
end
end
