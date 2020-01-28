%% ����rho�����е�NaN

function para = replot_rho(para, range)
rho = para.rho;
xi = linspace(range(3), range(4), 10000)';
rho = [xi, interp1(rho(:,1), rho(:,2), xi)];     % ����rho������һ�����޸�
rho(isnan(rho(:, 2)), 2) = 0;
para.rho = rho;
end