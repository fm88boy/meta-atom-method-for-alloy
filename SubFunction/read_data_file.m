%% 根据序号，读入势函数的数据

function [x, para, path, is_exist] = read_data_file(potential_id)
base_path = 'E:\Work\Matlab\PostWork';
x = []; para = [];

dir_path = sprintf('%s\\%d*', base_path, potential_id(1));
list = dir(dir_path);
base_path = sprintf('%s\\%s', base_path, list.name);

dir_path = sprintf('%s\\%d*', base_path, potential_id(2));
list = dir(dir_path);
path = sprintf('%s\\%s', base_path, list.name);

eam_path = sprintf('%s\\Work\\%d*.eam.fs', path, potential_id(3));
list = dir(eam_path);
if size(list,  1) > 0
    is_exist = 1;
else
    is_exist = 0;
end

filename = sprintf('%s\\Output_%d.mat', path, potential_id(3));
load(filename);
end