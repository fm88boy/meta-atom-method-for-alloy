%% 读取并初始化数据库

function conn = load_SQL()

dbname = 'potential';  % 调用数据库名称
username = 'fm88boy';
password = 'wp880112';

dataURL = sprintf('jdbc:mysql://localhost:3306/%s?useSSL=false', dbname);   % 把数据库放在这里面
conn = database('', username, password, 'com.mysql.jdbc.Driver', dataURL);
if ~isempty(conn.Message)
    error('数据库出错！')
end

my_exec(conn, 'set net_read_timeout=120');       % 提高数据库读取和写入的最大时间
my_exec(conn, 'set net_write_timeout=120');
my_exec(conn, 'set global max_allowed_packet=20*1024*1024');     % 提高数据库的最大命令程度为20Mb
% my_exec(conn, 'create database if not exists stock_index')  % 判断是否存在stock这个数据库，如果存在就使用这个数据库
% my_exec(conn, 'use stock_index');     % 使用这个数据库
end
