%% 消除rho函数中的NaN

function para = replot_rho(para, range)
rho = para.rho;
xi = linspace(range(3), range(4), 10000)';
rho = [xi, interp1(rho(:,1), rho(:,2), xi)];     % 对于rho函数做一定的修改
rho(isnan(rho(:, 2)), 2) = 0;
para.rho = rho;
end