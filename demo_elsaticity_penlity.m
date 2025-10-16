% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID = 'penlity_general_method.log';
logtype = 'wt';
% whether calculate compliance with respect to the same Emin
sameEmin = 5e-5;
same = false;
% parameters
volfrac = 0.3;          %体积比
nelx = 200;
nely = 100; 
Emin = [0.05,0.01,0.005,0.002];            %人造材料的相对杨氏模量
%bc = 'left_bdc_right_up_qin';
%bc = 'left_down_bdc_right_down_qin';%边界选
bc = 'left_up_down_bdc_left_down_qin';
%objectfunc = 'down_central';
%objectfunc = 'left_up';
objectfunc =  'right_down';
w = 4; %输出功权重
g = 0.000;
sd = 1;           %sd/nely 为卷积参数tau
continuation = 0; %是否使用预设的形状x
x = 1;            %预设形状x
V_constrain = 1;%1：不等式体积约束；0：等式体积约束
lambda = 1;
r      = 1000; %邻近因子
[y, loop, loop_k, c, x, energies, energies_k] = topthr_penlity_general(nelx, nely, lambda, r, volfrac, Emin(3), g, sd, objectfunc,bc, w,continuation, x, fileID,logtype,V_constrain);
%% 绘制目标函数收敛曲线对比
figure('Position', [100, 100, 800, 600]);  % 设置图形窗口大小

% 绘制不含线搜索的收敛曲线（红色+标记）
plot(1:loop, energies, 'r+', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName', '无线搜索');
hold on
% 绘制含线搜索的收敛曲线（绿色实线）
plot(1:loop_k, energies_k, 'g-', 'LineWidth', 2, ...
     'DisplayName', '含线搜索');

% 添加图例（清晰说明两种方法）
legend('不含线搜索次数', '含线搜索次数');

% 设置坐标轴标签
xlabel('迭代次数', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('目标函数值', 'FontSize', 14, 'FontWeight', 'bold');

% 添加标题
title('目标函数收敛曲线', ...
      'FontSize', 16, 'FontWeight', 'bold');

% 设置网格
grid on;
grid minor;

% 优化坐标轴
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;

% 添加文本注释（可选：解释收敛特性）
text(0.02, 0.98, ...
     {sprintf('最终收敛值: %.4f', energies(end))}, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', ...
     'FontSize', 10);
figure; imshow(1-x); 
