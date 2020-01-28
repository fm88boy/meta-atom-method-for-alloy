function [x0, para] = random_x_and_weight(x, para, rand_ratio)
x0 = x.*(1 + rand(size(x))*2*rand_ratio - rand_ratio);      % 对优化结果进行重新撒点
weight = para.weight;

name = fieldnames(weight);
for i = 1:size(name, 1)
    new_weight = getfield(weight, name{i}) * (1 + rand()*2*rand_ratio - rand_ratio);
    weight = setfield(weight, name{i}, new_weight);
end
weight = norm_struct(weight);       % 再将其归一化

para.weight = weight;
end