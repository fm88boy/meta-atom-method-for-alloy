%% ���ƺ������ΪEAM/FS��ʽ

function Write_EAM(F, rho, phi, name, element, tot_grid, filename, comment)      % ���ΪEAM�ļ�
fid = fopen(filename, 'w');
fprintf(fid, 'Produced by: Wang Peng, %s\nDate: %s\nContact information: wp_mech@zju.edu.cn\n', comment, datestr(now, 31));
fprintf(fid, '1  %s\n', char(name(1)));
fprintf(fid, '%d   %.15e   %d   %.15e   %.15e\n', tot_grid);
fprintf(fid, '%d   %.15e   %.15e   %s\n', element, char(name(2)));

temp = sprintf('%.15e  %.15e  %.15e  %.15e  %.15e\n', F(:, 2));
fprintf(fid, '%s', temp);

temp = sprintf('%.15e  %.15e  %.15e  %.15e  %.15e\n', rho(:, 2));
fprintf(fid, '%s', temp);

phi_mod = phi(:,1) .* phi(:, 2);    % ��FS�ļ����棬��r*phi�ĸ�ʽ
temp = sprintf('%.15e  %.15e  %.15e  %.15e  %.15e\n', phi_mod);
fprintf(fid, '%s', temp);

fclose(fid);
end
