%% 将一个结构体内部数据归一化

function data = norm_struct(data)
name = fieldnames(data);
total = 0;      % 记录了总数
for i = 1:size(name, 1)
    total = total + sum(getfield(data, name{i}));   % 获取值
end

for i = 1:size(name, 1)
    data = setfield(data, name{i}, getfield(data, name{i})/total);
end
end