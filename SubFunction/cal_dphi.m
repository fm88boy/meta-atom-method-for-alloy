% ���ݰ뾶���ڵ�ͽ״Σ������dphi����
function out = cal_dphi(radius, knots, order)
num_order = order(2) - order(1) + 1;    % ����һ���ڵ㣬�Ա����ĸ���
out = zeros(1, (size(knots, 2)-1)*num_order);         % ��ʼ��

for i = 1:size(out, 2)
    this_knots = ceil(i/num_order)+1;
    this_order = order(1) + mod(i-1, num_order);
    out(i) = heaviside(knots(this_knots)-radius)*(knots(this_knots)-radius).^(this_order-1)*(-this_order);
end
end
