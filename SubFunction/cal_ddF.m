% 根据rho计算出ddF的矩阵
function out = cal_ddF(rho)
out = [-0.25*rho^(-1.5), 2, 12*rho^2];
end
