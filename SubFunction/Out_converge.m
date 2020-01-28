% 输出参考文件
function Out_converge(out_txt, x, fval_new, step)
load Output_lowest para fval cons_iter
if fval_new < fval      % 如果这个值更小的话，那么就修改
    fval = fval_new;
    cons_iter = step;
%     para.base_data = x .* para.base_data;       % 修改模型的基础值
%     x = ones(size(x));
    save Output_lowest x para fval cons_iter;
end
lowest_step = cons_iter;
disp(out_txt);

fid = fopen('out.converge', 'a');
temp = sprintf('Date:\t%s', datestr(now, 31));
fprintf(fid, '%s\nStep:\t\t%d\tfval:\t\t%g\n', temp, step, fval_new);
fprintf(fid, 'Lowest_step:\t%d\tfval_lowest:\t%g\n', lowest_step, fval);
fprintf(fid, '%s\n', out_txt);
fclose(fid);
end
