% save the following output to this folder
save_to_folder = false;
% fileID for log file
logtype = 'wt';
descent_type = 'conjugate';
% whether calculate compliance with respect to the same Emin
sameEmin = 5e-5;
same = false;
% parameters
volfrac = 0.5;          %体积比
nelx = 400;
nely = 200; 
Emin = [0.05,0.01,0.005,0.001];            %人造材料的相对杨氏模量
bc = 'left_bdc';  %左边界固定
%bc = 'left_down_bdc_right_down_qin'; %右下角向左的推力
%bc = 'left_up_down_bdc_left_down_qin'; %左上角固定，左下角向右的推力
objectfunc = 'right_down_central'; %右边界中点向下的位移
%objectfunc = 'left_up'; %左上角向下的位移
%objectfunc =  'right_down'; %右下角向右的位移  
w = 1; %输出功权重
p = [1,0.8,0.6,0.4,0.2,-0.1,-0.3,-0.5,-0.7,-0.9,-1]; %GMIF参数
g = 0.000;
sd = 1;           %sd/nely 为卷积参数tau
continuation = 1; %是否使用预设的形状x
x_ini = zeros(nely,nelx);            %预设形状x
V_constrain = 1;%1：不等式体积约束；0：等式体积约束
lambda = 1;
xinitial = 1;
switch xinitial
    case 1 %中间一条1/5宽度的窄带
        len = floor(nely*volfrac); 
        lef = floor((nely-len)/2)+1;
        rig = lef + len -1;
        fixeddofs = [lef : rig]';
        ind = repmat(fixeddofs,1,nelx)+repmat((0:nelx-1)*nely, len, 1);
        Ind = reshape(ind, len*nelx,1);
    case 2
        len = floor(sqrt(nely*nelx*volfrac));
        lef = floor((nely-len)/2)+1;
        rig = lef + len -1;
        fixeddofs = [lef : rig]';
        ind = repmat(fixeddofs,1,len)+repmat(((lef-1):(rig-1))*nely, len, 1);
        Ind = reshape(ind, len^2,1);
    case 3
        len = floor(nely*volfrac/2); 
        lef = floor((nely-len)/2)+1;
        rig = lef + len -1;
        fixeddofs = [lef : rig]';
        ind = repmat(fixeddofs,1,nelx)+repmat((0:nelx-1)*nely, len, 1);
        Ind = union(reshape(ind, len*nelx,1), (((lef-1)*nely+1):(rig*nely))');
    case 4
        Ind = ones(nely,nelx);
end   
x_ini(Ind) = 1;
% if continuation == 1
%     figure; imshow(1-x);
% end
for k = 1:length(p)
r      = 10000000000; %邻近因子
fprintf('开始第%d次数值实验，p=%f.\n',k,p(k));
% fileID for log file
fileID = sprintf('%d-date.log', k);
[y, loop, loop_k, c, x, energies, energies_k] = topthr_penlity_general(nelx, nely, lambda, r, volfrac, Emin(3), p(k),g, sd, objectfunc,bc, w,continuation, x_ini, fileID,logtype,V_constrain,descent_type);
%% 绘制目标函数收敛曲线对比
%figure('Position', [100, 100, 800, 600]);  % 设置图形窗口大小

% 绘制不含线搜索的收敛曲线（红色+标记）
fig1 = figure('Name', sprintf('Figure 1 - Loop %d', k));
plot(1:loop, energies, 'r+', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName', 'Major iterations only');
hold on
% 绘制含线搜索的收敛曲线（绿色实线）
plot(1:loop_k, energies_k, 'g-', 'LineWidth', 2, ...
     'DisplayName', 'Including line search');

% 添加图例（清晰说明两种方法）
legend('Major iterations only', 'Including line search');

% 设置坐标轴标签
xlabel('Iterations', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Objective functional value', 'FontSize', 14, 'FontWeight', 'bold');

% 添加标题
title('Convergence history of the objective functional', ...
      'FontSize', 16, 'FontWeight', 'bold');

% 设置网格
grid on;
grid minor;

% 优化坐标轴
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;

% 添加文本注释（可选：解释收敛特性）
text(0.02, 0.98, ...
     {sprintf('Final value: %.4f', energies(end))}, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', ...
     'FontSize', 10);
%--------------------------------------------------
% 绘制不含线搜索的收敛曲线（红色+标记）
fig2 = figure('Name', sprintf('Figure 2 - Loop %d', k));
plot(1:loop, energies, 'r', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName','Major iterations only');
 % 设置坐标轴标签
xlabel('Iterations', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Objective functional value', 'FontSize', 14, 'FontWeight', 'bold');

% 添加标题
title('Convergence history of the objective functional', ...
      'FontSize', 16, 'FontWeight', 'bold');

% 设置网格
grid on;
grid minor;

% 优化坐标轴
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;

% 添加文本注释（可选：解释收敛特性）
text(0.02, 0.98, ...
     {sprintf('Final value: %.4f', energies(end))}, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', ...
     'FontSize', 10);
 
savefig(fig1, sprintf('%d-function_val_comp.fig', k));  % 保存FIG
print(fig1, sprintf('%d-function_val_comp.png', k), '-dpng', '-r300');  % 保存PNG，高分辨率

savefig(fig2, sprintf('%d-function_val.fig', k));
print(fig2, sprintf('%d-function_val.png', k), '-dpng', '-r300');

imwrite(1-x, sprintf('%d-shape.png', k), 'png');

% 关闭figure，避免窗口堆积（可选，如果想看可以注释掉）
close(fig1); close(fig2); 
end