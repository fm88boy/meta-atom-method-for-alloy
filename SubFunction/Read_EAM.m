%% 读入EAM/Alloy的合金势函数，可以读入多个元素，这里面把phi函数直接做了修正

function [F, rho, phi, N_tot, Element, tot_grid] = Read_EAM_Alloy(file_pair, num_row, r_span, rho_span, isplot)
% close all
% file_pair = 'E:\Temp.eam.alloy';
% r_span = [1.8, 6];
% rho_span = [0 4];
%% 读入EAM/Alloy文件，这里面考虑可能有多个原子

fid = fopen(file_pair);
for i = 1:3
    fgetl(fid);
end
temp = strtrim(fgetl(fid)); % 这里面包括本势函数包括元素的种类
N_tot = str2num(temp(1:2));    % 读入元素的个数，很明显最多不会超过99个元素
Element = cell(N_tot, 3);
temp = strtrim(temp(3:end));
for i = 1:N_tot-1
    for j = 1:size(temp, 2);
        if isspace(temp(j)) == 1
            break;
        end
    end
    Element{i, 1} = temp(1:j-1);
    temp = strtrim(temp(j:end));
end
Element{N_tot, 1} = strtrim(temp(1:end));      % 保存的是元素的种类

temp = strtrim(fgetl(fid));
tot_grid = sscanf(temp, '%d %f %d %f %f');  % 读入网格参数grid，依次为rho的网格数，rho的步长，phi的网格数，phi的步长，和cutoff

F = cell(N_tot, 1);
rho = cell(N_tot, 1);

for j = 1:N_tot
    temp = fgetl(fid);
    [element_single, ~, ~, next] = sscanf(temp, '%d %f %f');   % 读入元素的性质element，依次为元素序号，原子质量，晶格常数
    crystal = strtrim(temp(next:end));        % 读入构型晶格
    Element{j, 2} = element_single;
    Element{j, 3} = crystal;
    
    temp_F = zeros(tot_grid(1), 2);      % 读入F函数
    temp_F(:, 1) = 0 : tot_grid(2) : (tot_grid(2)*(tot_grid(1)-1));
    for i = 1 : ceil(tot_grid(1)/num_row)
        temp = lower(fgetl(fid));
        t_data = str2num(temp);
        t_start = (i-1)*num_row+1;
        t_end = i*num_row;
        temp_F(t_start:t_end, 2) = t_data;
    end
    F{j} = temp_F;
    
    temp_rho = zeros(tot_grid(3), 2);    % 读入rho函数
    temp_rho(:, 1) = 0 : tot_grid(4) : (tot_grid(4)*(tot_grid(3)-1));
    for i = 1 : ceil(tot_grid(3)/num_row)
        temp = lower(fgetl(fid));
        t_data = str2num(temp);
        t_start = (i-1)*num_row+1;
        t_end = i*num_row;
        temp_rho(t_start:t_end, 2) = t_data;
    end
    rho{j} = temp_rho;
end

phi = cell(N_tot*(N_tot+1)/2, 1);   % 考虑所有可能的组合
for j = 1:size(phi, 1)
    temp_phi = zeros(tot_grid(3), 2);    % 读入phi函数
    temp_phi(:, 1) = 0:tot_grid(4):(tot_grid(4)*(tot_grid(3)-1));
    for i = 1 : tot_grid(3)/num_row
%         fprintf('%d\n', i)
        temp = lower(fgetl(fid));
        t_data = str2num(temp);
        t_start = (i-1)*num_row+1;
        t_end = i*num_row;
        temp_phi(t_start:t_end, 2) = t_data;
    end
    temp_phi(:,2) = temp_phi(:,2)./temp_phi(:,1);
    phi{j} = temp_phi;
end
fclose(fid);

%% 分析数据
temp = sprintf('Number of elements: %d', N_tot);
disp(temp)

temp_str = '';
for i = 1:N_tot
    temp_str = sprintf('%s %s', temp_str, Element{i, 1});
end
temp = sprintf('Lists of elements: %s', temp_str);
disp(temp)

if isplot == 1
    str_F = cell(1, N_tot);
    str_rho = cell(1, N_tot);
    close all
    for i = 1:N_tot
        figure(1)
        hold on
        F_single = F{i};
        plot(F_single(:, 1), F_single(:, 2))
        str_F{i} = sprintf('F: %s', Element{i, 1});
        hold off
        
        figure(2)
        hold on
        rho_single = rho{i};
        plot(rho_single(:, 1), rho_single(:, 2))
        str_rho{i} = sprintf('rho: %s', Element{i, 1});
        hold off
    end
    figure(1)
    title('F')
    legend(str_F)
    xlim(rho_span)
    
    figure(2)
    title('rho')
    legend(str_rho)
    xlim(r_span)

    index_1 = 1;
    index_2 = 1;
    str_phi = cell(1, size(phi, 1));
    for i = 1:size(phi, 1)
        figure(3)
        hold on
        str_phi{i} = sprintf('phi: %s-%s', Element{index_1, 1}, Element{index_2, 1});
        phi_single = phi{i};
        plot(phi_single(:, 1), phi_single(:, 2))
        hold off
        
        if index_2 < index_1
            index_2 = index_2 + 1;
        else
            index_1 = index_1 + 1;
            index_2 = 1;
        end
    end
    
    figure(3)
    title('phi')
    legend(str_phi)
    xlim(r_span)    
end
% save EAM_alloy F rho phi N_tot Element tot_grid
end
