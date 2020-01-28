% ���ݰ뾶���ڵ�ͽ״Σ������phi����
function out = cal_phi(radius, knots, order)
% radius��ԭ�Ӱ뾶��knots���ڵ㣻order������ʽ�״�
num_order = order(2) - order(1) + 1;    % ����һ���ڵ㣬�Ա����ĸ���
out = zeros(1, (size(knots, 2)-1)*num_order);         % ��һ��phi��knots������ȷ�����޵�

for i = 1:size(out, 2)
    this_knots = ceil(i/num_order)+1;   % ��Ӧ�Ľڵ�
    this_order = order(1) + mod(i-1, num_order);    % ��Ӧ�Ľ״�
    out(i) = heaviside(knots(this_knots)-radius)*(knots(this_knots)-radius).^this_order;
end
end
