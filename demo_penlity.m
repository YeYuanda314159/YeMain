% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID = 'D:\Study\thresholding dynamics\Numerical Results2\24  -data.log';
% whether calculate compliance with respect to the same Emin
% parameters
nelx = 600; 
nely = 600;
volfrac = 0.2; %体积占比
lambda = 0.01; %正则参数
p      = 0.5; %隐式罚参数
r      = 1000; %邻近因子
g = 0.00001; %周长罚参数
descent_type = 'conjugate';
sd = 1;
bc = 'left_Dirichlet'; %左端1/5Dirichlet边界
%bc = 'all_Dirichlet'; %完全Dirichlet边界
%bc = 'topleft_Dirichlet';
%bc = 'allleft_Dirichlet';
continuation = 1; %是否使用预设初始形状，0：默认均匀初始值; 1：使用给定初始值
V_constrain = 0; %0:等式体积约束，1：不等式体积约束
xinitial = 2;
x = zeros(nely,nelx);
switch xinitial
    case 1 %中间一条1/5宽度的窄带
        len = floor(nely*volfrac); 
        lef = (nely-len)/2+1;
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
        lef = (nely-len)/2+1;
        rig = lef + len -1;
        fixeddofs = [lef : rig]';
        ind = repmat(fixeddofs,1,nelx)+repmat((0:nelx-1)*nely, len, 1);
        Ind = union(reshape(ind, len*nelx,1), (((lef-1)*nely+1):(rig*nely))');
end   
if continuation == 1
    x(Ind) = 1;
    figure; imshow(1-x);
end
[y,loop,loop_k,c,x,energies,energies_k]=topthr_penlity(nelx,nely,volfrac,lambda,p,r,g,sd,bc,continuation,x,fileID,V_constrain,descent_type);
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
title('目标函数收敛曲线对比', ...
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
%---------------------------------
figure('Position', [100, 100, 800, 600]);  % 设置图形窗口大小

% 绘制不含线搜索的收敛曲线（红色+标记）
plot(1:loop, energies, 'r', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName', '无线搜索');
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