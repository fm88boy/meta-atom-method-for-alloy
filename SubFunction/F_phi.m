%% ���ݴ��뺯���Ĳ�������ԭ��F������phi����

function [F, phi] = F_phi(range, data, split, knots, order, N)
% range��F������phi������ȡֵ��Χ��data��F������phi�����Ĳ���ֵ
% split��F������phi�����ķָ�㣻knots��phi�����Ľڵ㣬order��phi�����Ľ״�
% N��F������phi�������Ա������ݵ�

F_data = data(1:split-1);   % �����ݷֿ�
phi_data = data(split:end);

xi = linspace(range(1), range(2), N)';
yi = F_data(1)*xi.^(0.5) + F_data(2)*xi.^2 + F_data(3)*xi.^4;
F = [xi, yi];
F(:,2) = F(:,2) + heaviside(F(:,1)-60) * 0.2 .* (F(:,1)-60).^2;    % ���������Ҹ���һ��������

% ����[0, knots(1)]������Ϊa*x^2+bx+c
y1 = dot(cal_ddphi(knots(1), knots, order), phi_data)/2;
y2 = dot(cal_dphi(knots(1), knots, order), phi_data) - 2*y1*knots(1);
y3 = dot(cal_phi(knots(1), knots, order), phi_data) - y1*(knots(1))^2 - y2*knots(1);
% ��������Ϊa,b,c��ȷ����knots(1)��������������
coeff = [y1, y2, y3];
inner_phi = coeff(1)*0.5.^2 + coeff(2)*0.5 + coeff(3);

xi = linspace(range(3), range(4), N)';
sum_temp = zeros(size(xi));         % ��ʼ��
num_order = order(2) - order(1) + 1;    % ����һ���ڵ㣬�Ա����ĸ���
for i = 2:size(knots, 2)
    for j = order(1) : order(2)
        % ����Ӧ����ȡ���������Ҽ��뵽phi������
        sum_temp = sum_temp + heaviside(knots(i)-xi).*(knots(i)-xi).^j.*phi_data((i-2)*num_order+j-order(1)+1);
    end
end
yi = heaviside(0.5-xi).*inner_phi + heaviside(xi-0.5).*heaviside(knots(1)-xi).*(coeff(1)*xi.^2+coeff(2)*xi+coeff(3)) ...
    + heaviside(xi-knots(1)).*sum_temp;
phi = [xi, yi];
end
