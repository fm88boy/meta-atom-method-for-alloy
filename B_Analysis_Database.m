%% �����ݿ��и���Ҫ���ȡ�ƺ���

function B_Analysis_Database()
initialize();
table_name = 'HEA_Fe50Mn30Cr10Co10';
element = {'HEA', 'fcc', 26, 55.4968};        % Ԫ����ź�ԭ������
max_num = 50;                                 % ��������
tolerance.least_fval = 10;                     % ��Ϊ�㣬��Ԣ���Ŵӿɲ����ƺ�����ѡȡ���в���С��N���ƺ���
tolerance.id_range = [20066000, 1e9];         % ɨ��ķ�Χ
tolerance.lattice = 0.1;
tolerance.c_a_ratio = 0.5;
tolerance.elastic = 0.2;
tolerance.sub_energy = 0.1;
tolerance.sf_energy = 0.1;
tolerance.usf_energy = 0.5;
tolerance.surface_energy = 0.5;
tolerance.vac_energy = 1;

tolerance.fval_min = 10;            % ��Ͳв�Ҫ��

tolerance.fcc_strain_100 = 0;          % 0������ֻҪ�����е�·���������0��δ��Ҳ�����������ɸѡ����
tolerance.fcc_strain_110 = 0;          % 0������ֻҪ�����е�·���������0��δ��Ҳ�����������ɸѡ����
tolerance.fcc_strain_111 = 0;          % 0������ֻҪ�����е�·���������0��δ��Ҳ�����������ɸѡ����
tolerance.fcc_specific_energy = 0;     % 0������ֻҪ�����е�·���������0��δ��Ҳ�����������ɸѡ����

is_delete = 0;                      % �Ƿ�ɾ���ƺ���

% cons_select = [19, 18];           % ѡ��Լ����������Դ�����Ϊ-1����ѡȡ���һ������
cons_select = -1;                   % ѡ��Լ����������Դ�����Ϊ-1����ѡȡ���һ������

%% ��ʼ��ȡ���ݿ�
conn = load_SQL();
% my_exec(conn, 'create table if not exists HEA_Fe50Mn30Cr10Co10 (id INT unsigned not null primary key, group_id INT not null, task_id INT not null, potential_id INT not null, fval FLOAT not null, is_test INT not null, a_fcc FLOAT, C11_fcc FLOAT, C12_fcc FLOAT, C44_fcc FLOAT, esub_fcc FLOAT, esf_fcc FLOAT, eusf_fcc FLOAT, e111_fcc FLOAT, e110_fcc FLOAT, e100_fcc FLOAT, evac_fcc FLOAT, a_hcp FLOAT, c_hcp FLOAT, esub_hcp FLOAT, Bain_1_strain FLOAT, Bain_1_energy FLOAT, Bain_2_strain FLOAT, Bain_2_energy FLOAT, Bain_3_strain FLOAT, Bain_3_energy FLOAT, Bain_4_strain FLOAT, Bain_4_energy FLOAT, phi_diff_energy FLOAT, max_phi FLOAT)')
% my_exec(conn, 'create table if not exists HEA_Fe50Mn30Cr10Co10_cons (id INT unsigned not null primary key, group_id INT not null, task_id INT not null, a_fcc FLOAT, C11_fcc FLOAT, C12_fcc FLOAT, C44_fcc FLOAT, esub_fcc FLOAT, esf_fcc FLOAT, eusf_fcc FLOAT, e111_fcc FLOAT, e110_fcc FLOAT, e100_fcc FLOAT, evac_fcc FLOAT, a_hcp FLOAT, c_hcp FLOAT, esub_hcp FLOAT, Bain_1_strain FLOAT, Bain_1_energy_low FLOAT, Bain_1_energy_high FLOAT, Bain_2_strain FLOAT, Bain_2_energy_low FLOAT, Bain_2_energy_high FLOAT, Bain_3_strain FLOAT, Bain_3_energy_low FLOAT, Bain_3_energy_high FLOAT, Bain_4_strain FLOAT, Bain_4_energy_low FLOAT, Bain_4_energy_high FLOAT, phi_diff_energy FLOAT, max_phi FLOAT)')

[own_index, ~] = read_potential(conn, cons_select, table_name, tolerance);
fprintf('�ƺ���������%d\n', size(own_index, 1))

if is_delete == 1           % �����ȷ����ɾ���ƺ���
    for i = 1:size(own_index, 1)
        fprintf('ID: %d-%d-%d, residual error: %.4f\n', own_index(i, 2:5));
        answer = input('�����ɾ����������y���س����������ɾ������ֱ�ӻس���', 's');
        if strcmp(answer, 'y')
            id = own_index(i, 2)*1e6 + own_index(i, 3)*1e3 + own_index(i, 4);
            command = sprintf('update %s set is_test = -1 where id = %d', table_name, id);
            my_exec(conn, command);
        end
    end
    close(conn);
    return;
end

if size(own_index, 1) <= max_num
    for i = 1:size(own_index, 1)
        [x, para, path, is_exist] = read_data_file(own_index(i, 2:4));
        write_para.path = path;
        write_para.comment = sprintf('ID: %d-%d-%d, residual error: %.4f', own_index(i, 2:5));
        write_para.element = element;
        write_para.write_index = own_index(i, 4);
        disp(write_para.comment);
        
        if is_exist == 0
            [out_txt, ~] = Post_analysis(x, para, 1, write_para);
            disp(out_txt);
        end
    end
elseif size(own_index, 1) > max_num
    disp('����Լ������');
end
close(conn);

end

function initialize()
sub_path = sprintf('%s/SubFunction', pwd);
addpath(sub_path);
end

%% �����ƺ���
function [own_index, potential_cons] = read_potential(conn, cons_select, table_name, tolerance)
table_name_cons = strcat(table_name, '_cons');          % ���ȶ���Լ����
temp{1} = 'a_fcc, C11_fcc, C12_fcc, C44_fcc, esub_fcc, esf_fcc, eusf_fcc, e111_fcc, e110_fcc, e100_fcc, evac_fcc';
temp{2} = 'a_hcp, c_hcp, esub_hcp, a_bcc, esub_bcc';
label_name = sprintf('%s, ', temp{:});
label_name = label_name(1:end-2);

if cons_select(1) == -1
    command = sprintf('select %s from %s order by id desc limit 1', label_name, table_name_cons);
    data_cons = my_exec(conn, command);
else
    command = sprintf('select %s from %s where group_id = %d and task_id = %d', label_name, table_name_cons, cons_select(1), cons_select(2));
    data_cons = my_exec(conn, command);
end
potential_cons = cell2mat(data_cons);

own_index = read_potential_inner(conn, potential_cons, table_name, tolerance);
end

function own_index = read_potential_inner(conn, potential_cons, table_name, tolerance)
if tolerance.least_fval == 0
    % һ��14������
    a_fcc = potential_cons(1); C11_fcc = potential_cons(2); C12_fcc = potential_cons(3); C44_fcc = potential_cons(4);
    esub_fcc = potential_cons(5); esf_fcc = potential_cons(6); eusf_fcc = potential_cons(7); e111_fcc = potential_cons(8);
    e110_fcc = potential_cons(9); e100_fcc = potential_cons(10); evac_fcc = potential_cons(11); a_hcp = potential_cons(12);
    c_hcp = potential_cons(13); esub_hcp = potential_cons(14); a_bcc = potential_cons(15); esub_bcc = potential_cons(16);
    
    range.a_fcc = sort(a_fcc * (1 + [-1, 1]*tolerance.lattice));
    range.C11_fcc = sort(C11_fcc * (1 + [-1, 1]*tolerance.elastic));
    range.C12_fcc = sort(C12_fcc * (1 + [-1, 1]*tolerance.elastic));
    range.C44_fcc = sort(C44_fcc * (1 + [-1, 1]*tolerance.elastic));
    range.esub_fcc = sort(esub_fcc * (1 + [-1, 1]*tolerance.sub_energy));
    range.esf_fcc = sort(esf_fcc * (1 + [-1, 1]*tolerance.sf_energy));
    range.eusf_fcc = sort(eusf_fcc * (1 + [-1, 1]*tolerance.usf_energy));
    range.e111_fcc = sort(e111_fcc * (1 + [-1, 1]*tolerance.surface_energy));
    range.e110_fcc = sort(e110_fcc * (1 + [-1, 1]*tolerance.surface_energy));
    range.e100_fcc = sort(e100_fcc * (1 + [-1, 1]*tolerance.surface_energy));
    range.evac_fcc = sort(evac_fcc * (1 + [-1, 1]*tolerance.vac_energy));
    % HCP��صĲ���
    range.a_hcp = sort(a_hcp * (1 + [-1, 1]*tolerance.lattice));
    range.c_hcp = sort(c_hcp * (1 + [-1, 1]*tolerance.c_a_ratio));
    range.esub_hcp = sort(esub_hcp * (1 + [-1, 1]*tolerance.sub_energy));
    % BCC��صĲ���
    range.a_bcc = sort(a_bcc * (1 + [-1, 1]*tolerance.lattice));
    range.esub_bcc = sort(esub_bcc * (1 + [-1, 1]*tolerance.sub_energy));
    
    temp{1} = sprintf('a_fcc between %.6f and %.6f and C11_fcc between %.6f and %.6f and C12_fcc between %.6f and %.6f and C44_fcc between %.6f and %.6f', ...
        range.a_fcc, range.C11_fcc, range.C12_fcc, range.C44_fcc);
    temp{2} = sprintf('esub_fcc between %.6f and %.6f and esf_fcc between %.6f and %.6f and eusf_fcc between %.6f and %.6f and e111_fcc between %.6f and %.6f', ...
        range.esub_fcc, range.esf_fcc, range.eusf_fcc, range.e111_fcc);
    temp{3} = sprintf('e110_fcc between %.6f and %.6f and e100_fcc between %.6f and %.6f and evac_fcc between %.6f and %.6f and a_hcp between %.6f and %.6f', ...
        range.e110_fcc, range.e100_fcc, range.evac_fcc, range.a_hcp);
    temp{4} = sprintf('c_hcp between %.6f and %.6f and esub_hcp between %.6f and %.6f and a_bcc between %.6f and %.6f and esub_bcc between %.6f and %.6f', ...
        range.c_hcp, range.esub_hcp, range.a_bcc, range.esub_bcc);
    constraint = sprintf('%s and ', temp{:});
    constraint = constraint(1:end-5);
    
    command = sprintf('select * from %s where fval < %.6f and is_test >=0 and %s', table_name, tolerance.fval_min, constraint);
    data = my_exec(conn, command);
    if ~strcmp(data{1}, 'No Data')
        own_index = cell2mat(data);
    else
        own_index = [];
    end
else
    command = sprintf('select * from %s where is_test >=0 and id between %d and %d order by fval limit %d', ...
        table_name, tolerance.id_range(1), tolerance.id_range(2), tolerance.least_fval);
    data = my_exec(conn, command);
    if ~strcmp(data{1}, 'No Data')
        own_index = cell2mat(data);
    else
        own_index = [];
    end    
end
end
