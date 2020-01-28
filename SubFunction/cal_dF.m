% 根据rho计算出dF的矩阵
function out = cal_dF(rho)
out = [0.5*rho^(-0.5), 2*rho, 4*rho^3];
end

