%% ��ȡ����ʼ�����ݿ�

function conn = load_SQL()

dbname = 'potential';  % �������ݿ�����
username = 'fm88boy';
password = 'wp880112';

dataURL = sprintf('jdbc:mysql://localhost:3306/%s?useSSL=false', dbname);   % �����ݿ����������
conn = database('', username, password, 'com.mysql.jdbc.Driver', dataURL);
if ~isempty(conn.Message)
    error('���ݿ����')
end

my_exec(conn, 'set net_read_timeout=120');       % ������ݿ��ȡ��д������ʱ��
my_exec(conn, 'set net_write_timeout=120');
my_exec(conn, 'set global max_allowed_packet=20*1024*1024');     % ������ݿ���������̶�Ϊ20Mb
% my_exec(conn, 'create database if not exists stock_index')  % �ж��Ƿ����stock������ݿ⣬������ھ�ʹ��������ݿ�
% my_exec(conn, 'use stock_index');     % ʹ��������ݿ�
end
