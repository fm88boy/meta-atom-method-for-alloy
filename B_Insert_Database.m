%% 读取势函数的结果，进行合理分析后，保存到数据库中

function B_Insert_Database()
initialize();
Path = 'E:\Work\Matlab\PostWork\20_Fe50Mn30Cr10Co10_3\';
search = '*D1222*';     % 搜索的关键字
table_name = 'HEA_Fe50Mn30Cr10Co10';
lowest_fval = 1;

%% 开始读取数据库
conn = load_SQL();
recycle('on')   % 保存到回收站

% my_exec(conn, 'create table if not exists HEA_Fe50Mn30Cr10Co10 (id INT unsigned not null primary key, group_id INT not null, task_id INT not null, potential_id INT not null, fval FLOAT not null, a_fcc FLOAT, C11_fcc FLOAT, C12_fcc FLOAT, C44_fcc FLOAT, esub_fcc FLOAT, esf_fcc FLOAT, eusf_fcc FLOAT, e111_fcc FLOAT, e110_fcc FLOAT, e100_fcc FLOAT, evac_fcc FLOAT, a_hcp FLOAT, c_hcp FLOAT, esub_hcp FLOAT, Bain_1_strain FLOAT, Bain_1_energy FLOAT, Bain_2_strain FLOAT, Bain_2_energy FLOAT, Bain_3_strain FLOAT, Bain_3_energy FLOAT, Bain_4_strain FLOAT, Bain_4_energy FLOAT, phi_diff_energy FLOAT, max_phi FLOAT)')
% my_exec(conn, 'create table if not exists HEA_Fe50Mn30Cr10Co10_cons (id INT unsigned not null primary key, group_id INT not null, task_id INT not null, a_fcc FLOAT, C11_fcc FLOAT, C12_fcc FLOAT, C44_fcc FLOAT, esub_fcc FLOAT, esf_fcc FLOAT, eusf_fcc FLOAT, e111_fcc FLOAT, e110_fcc FLOAT, e100_fcc FLOAT, evac_fcc FLOAT, a_hcp FLOAT, c_hcp FLOAT, esub_hcp FLOAT, Bain_1_strain FLOAT, Bain_1_energy_low FLOAT, Bain_1_energy_high FLOAT, Bain_2_strain FLOAT, Bain_2_energy_low FLOAT, Bain_2_energy_high FLOAT, Bain_3_strain FLOAT, Bain_3_energy_low FLOAT, Bain_3_energy_high FLOAT, Bain_4_strain FLOAT, Bain_4_energy_low FLOAT, Bain_4_energy_high FLOAT, phi_diff_energy FLOAT, max_phi FLOAT)')

temp = regexp(Path, '\\(\d*)_', 'tokens');
if ~isempty(temp)
    group_id = str2double(temp{1});     % group_id：分类
end
table_name_cons = strcat(table_name, '_cons');

list = dir(strcat(Path, search));       % 外部路径
for i = 1:size(list, 1)
    temp = regexp(list(i).name, '(\d*)_D', 'tokens');
    if ~isempty(temp)
        task_id = str2double(temp{1});  % task_id: 任务
    end

    inner_Path = strcat(Path, list(i).name, '\');
    inner_list = dir(strcat(inner_Path, 'Output_*.mat'));           % 扫描目录内的势函数文件
    for j = 1:size(inner_list, 1)
        filename = inner_list(j).name;
        temp = regexp(filename, 'Output_(\d*).mat', 'tokens');
        full_pathname = strcat(inner_Path, filename);
        if ~isempty(temp)
            potential_id = str2double(temp{1});
            load(full_pathname);
        else
            continue;           % 如果不是正常的数据，则跳过
        end
        [~, para_out] = Post_analysis(x, para, 0, []);
        if potential_id == 1        % 如果是第一个势函数，则保存约束量
            insert_cons(conn, group_id, task_id, potential_id, table_name_cons, para_out);
        end
        insert_data(conn, group_id, task_id, potential_id, table_name, para_out);
        fprintf('已经读入%d-%d-%d\n', group_id, task_id, potential_id);
        if para_out.out > lowest_fval       % 如果残差值过大，则删除该文件
            delete(full_pathname)
            disp('删除文件。')
        end
    end
end
close(conn);
end

function initialize()
sub_path = sprintf('%s/SubFunction', pwd);
addpath(sub_path);
end

function insert_cons(conn, group_id, task_id, potential_id, table_name, para)
id = group_id*1e6 + task_id*1e3 + potential_id;
a_fcc = para.cons_data.lattice(1);
a_hcp = para.cons_data.lattice(2);
c_hcp = para.cons_data.lattice(3);
a_bcc = para.cons_data.lattice(4);
C11_fcc = para.cons_data.fcc_elastic(1);
C12_fcc = para.cons_data.fcc_elastic(2);
C44_fcc = para.cons_data.fcc_elastic(3);
esub_fcc = para.cons_data.sub_energy(1);
esub_hcp = para.cons_data.sub_energy(2);
esub_bcc = para.cons_data.sub_energy(3);
esf_fcc = para.cons_data.sf_energy(1);
eusf_fcc = para.cons_data.usf_energy(1);
e111_fcc = para.cons_data.fcc_surface_energy(1);
e110_fcc = para.cons_data.fcc_surface_energy(2);
e100_fcc = para.cons_data.fcc_surface_energy(3);
evac_fcc = para.cons_data.vac_energy(1);

fcc_strain_100 = para.cons_data.fcc_strain_100;
fcc_strain_110 = para.cons_data.fcc_strain_110;
fcc_strain_111 = para.cons_data.fcc_strain_111;
phi_diff_energy = para.cons_data.phi_well(4);
max_phi = para.cons_data.phi_well(6);
fcc_specific_energy = para.cons_data.fcc_specific_energy;

temp_1{1} = 'id, group_id, task_id';
temp_1{2} = 'a_fcc, C11_fcc, C12_fcc, C44_fcc, esub_fcc, esf_fcc, eusf_fcc, e111_fcc, e110_fcc, e100_fcc, evac_fcc';
temp_1{3} = 'a_hcp, c_hcp, esub_hcp, a_bcc, esub_bcc';
label_name = sprintf('%s, ', temp_1{:});
label_name = label_name(1:end-2);

temp_2{1} = sprintf('%d, %d, %d', id, group_id, task_id);
temp_2{2} = sprintf('%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f', ...
    a_fcc, C11_fcc, C12_fcc, C44_fcc, esub_fcc, esf_fcc, eusf_fcc, e111_fcc, e110_fcc, e100_fcc, evac_fcc);
temp_2{3} = sprintf('%.6f, %.6f, %.6f, %.6f, %.6f', ...
    a_hcp, c_hcp, esub_hcp, a_bcc, esub_bcc);
label_value = sprintf('%s, ', temp_2{:});
label_value = label_value(1:end-2);

command = sprintf('replace into %s (%s) values (%s)', table_name, label_name, label_value);
my_exec(conn, command);

table_name = strcat(table_name, '_addition');       % 附加条件
temp_3{1} = 'id, group_id, task_id';
temp_3{2} = 'fcc_100_strain, fcc_110_strain, fcc_111_strain, phi_diff_energy, max_phi, fcc_specific_energy';
str_fcc_100_strain = encode_matrix_to_string(fcc_strain_100);
str_fcc_110_strain = encode_matrix_to_string(fcc_strain_110);
str_fcc_111_strain = encode_matrix_to_string(fcc_strain_111);
str_fcc_specific_energy = encode_matrix_to_string(fcc_specific_energy);
label_name = sprintf('%s, ', temp_3{:});
label_name = label_name(1:end-2);

temp_4{1} = sprintf('%d, %d, %d', id, group_id, task_id);
temp_4{2} = sprintf('''%s'', ''%s'', ''%s'', %.6f, %.6f, ''%s''', ...
    str_fcc_100_strain, str_fcc_110_strain, str_fcc_111_strain, phi_diff_energy, max_phi, str_fcc_specific_energy);
label_value = sprintf('%s, ', temp_4{:});
label_value = label_value(1:end-2);
command = sprintf('replace into %s (%s) values (%s)', table_name, label_name, label_value);
my_exec(conn, command);
end

function insert_data(conn, group_id, task_id, potential_id, table_name, para)
id = group_id*1e6 + task_id*1e3 + potential_id;
fval = para.out;        % 基础残差
is_test = 0;            % 0代表着没有经过测试
a_fcc = para.cal_data.lattice(1);
a_hcp = para.cal_data.lattice(2);
c_hcp = para.cal_data.lattice(3);
a_bcc = para.cal_data.lattice(4);
C11_fcc = para.cal_data.fcc_elastic(1);
C12_fcc = para.cal_data.fcc_elastic(2);
C44_fcc = para.cal_data.fcc_elastic(3);
esub_fcc = para.cal_data.sub_energy(1);
esub_hcp = para.cal_data.sub_energy(2);
esub_bcc = para.cal_data.sub_energy(3);
esf_fcc = para.cal_data.sf_energy(1);
eusf_fcc = para.cal_data.usf_energy(1);
e111_fcc = para.cal_data.fcc_surface_energy(1);
e110_fcc = para.cal_data.fcc_surface_energy(2);
e100_fcc = para.cal_data.fcc_surface_energy(3);
evac_fcc = para.cal_data.vac_energy(1);

fcc_strain_100 = para.cal_data.fcc_strain_100;
fcc_strain_110 = para.cal_data.fcc_strain_110;
fcc_strain_111 = para.cal_data.fcc_strain_111;
phi_diff_energy = para.cal_data.phi_well(1);
max_phi = para.cal_data.phi_well(2);
fcc_specific_energy = para.cal_data.fcc_specific_energy;

temp_1{1} = 'id, group_id, task_id, potential_id, fval, is_test';
temp_1{2} = 'a_fcc, C11_fcc, C12_fcc, C44_fcc, esub_fcc, esf_fcc, eusf_fcc, e111_fcc, e110_fcc, e100_fcc, evac_fcc';
temp_1{3} = 'a_hcp, c_hcp, esub_hcp, a_bcc, esub_bcc';
label_name = sprintf('%s, ', temp_1{:});
label_name = label_name(1:end-2);

temp_2{1} = sprintf('%d, %d, %d, %d, %.6f, %d', id, group_id, task_id, potential_id, fval, is_test);
temp_2{2} = sprintf('%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f', ...
    a_fcc, C11_fcc, C12_fcc, C44_fcc, esub_fcc, esf_fcc, eusf_fcc, e111_fcc, e110_fcc, e100_fcc, evac_fcc);
temp_2{3} = sprintf('%.6f, %.6f, %.6f, %.6f, %.6f', a_hcp, c_hcp, esub_hcp, a_bcc, esub_bcc);
label_value = sprintf('%s, ', temp_2{:});
label_value = label_value(1:end-2);

command = sprintf('replace into %s (%s) values (%s)', table_name, label_name, label_value);
my_exec(conn, command);

table_name = strcat(table_name, '_addition');       % 附加条件
temp_3{1} = 'id, group_id, task_id, potential_id';
temp_3{2} = 'fcc_100_strain, fcc_110_strain, fcc_111_strain, phi_diff_energy, max_phi, fcc_specific_energy';
str_fcc_100_strain = encode_matrix_to_string(fcc_strain_100);
str_fcc_110_strain = encode_matrix_to_string(fcc_strain_110);
str_fcc_111_strain = encode_matrix_to_string(fcc_strain_111);
str_fcc_specific_energy = encode_matrix_to_string(fcc_specific_energy);
label_name = sprintf('%s, ', temp_3{:});
label_name = label_name(1:end-2);

temp_4{1} = sprintf('%d, %d, %d, %d', id, group_id, task_id, potential_id);
temp_4{2} = sprintf('''%s'', ''%s'', ''%s'', %.6f, %.6f, ''%s''', ...
    str_fcc_100_strain, str_fcc_110_strain, str_fcc_111_strain, phi_diff_energy, max_phi, str_fcc_specific_energy);
label_value = sprintf('%s, ', temp_4{:});
label_value = label_value(1:end-2);
command = sprintf('replace into %s (%s) values (%s)', table_name, label_name, label_value);
my_exec(conn, command);
end
