%% ��һ���ṹ���ڲ����ݹ�һ��

function data = norm_struct(data)
name = fieldnames(data);
total = 0;      % ��¼������
for i = 1:size(name, 1)
    total = total + sum(getfield(data, name{i}));   % ��ȡֵ
end

for i = 1:size(name, 1)
    data = setfield(data, name{i}, getfield(data, name{i})/total);
end
end