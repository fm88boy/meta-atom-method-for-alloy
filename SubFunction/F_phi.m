%% 根据传入函数的参数，复原出F函数和phi函数

function [F, phi] = F_phi(range, data, split, knots, order, N)
% range：F函数和phi函数的取值范围；data：F函数和phi函数的参数值
% split：F函数和phi函数的分割点；knots：phi函数的节点，order：phi函数的阶次
% N：F函数和phi函数的自变量数据点

F_data = data(1:split-1);   % 把数据分开
phi_data = data(split:end);

xi = linspace(range(1), range(2), N)';
yi = F_data(1)*xi.^(0.5) + F_data(2)*xi.^2 + F_data(3)*xi.^4;
F = [xi, yi];
F(:,2) = F(:,2) + heaviside(F(:,1)-60) * 0.2 .* (F(:,1)-60).^2;    % 在这里面我给了一个修正项

% 对于[0, knots(1)]，函数为a*x^2+bx+c
y1 = dot(cal_ddphi(knots(1), knots, order), phi_data)/2;
y2 = dot(cal_dphi(knots(1), knots, order), phi_data) - 2*y1*knots(1);
y3 = dot(cal_phi(knots(1), knots, order), phi_data) - y1*(knots(1))^2 - y2*knots(1);
% 这三个量为a,b,c，确保在knots(1)处函数二阶连续
coeff = [y1, y2, y3];
inner_phi = coeff(1)*0.5.^2 + coeff(2)*0.5 + coeff(3);

xi = linspace(range(3), range(4), N)';
sum_temp = zeros(size(xi));         % 初始化
num_order = order(2) - order(1) + 1;    % 对于一个节点，自变量的个数
for i = 2:size(knots, 2)
    for j = order(1) : order(2)
        % 将对应的量取出来，并且加入到phi函数中
        sum_temp = sum_temp + heaviside(knots(i)-xi).*(knots(i)-xi).^j.*phi_data((i-2)*num_order+j-order(1)+1);
    end
end
yi = heaviside(0.5-xi).*inner_phi + heaviside(xi-0.5).*heaviside(knots(1)-xi).*(coeff(1)*xi.^2+coeff(2)*xi+coeff(3)) ...
    + heaviside(xi-knots(1)).*sum_temp;
phi = [xi, yi];
end
