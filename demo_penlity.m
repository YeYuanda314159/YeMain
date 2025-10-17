% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID = 'D:\Study\thresholding dynamics\ICTM-split-penlity\example1.log';
% whether calculate compliance with respect to the same Emin
% parameters
nelx = 200; 
nely = 200;
volfrac = 0.2; %体积占比
lambda = 10000; %正则参数
r      = 100; %邻近因子
g = 0.00001; %周长罚参数
sd = 1;
bc = 'left_Dirichlet'; %左端1/5Dirichlet边界
%bc = 'all_Dirichlet'; %完全Dirichlet边界
%bc = 'topleft_Dirichlet';
%bc = 'allleft_Dirichlet';
continuation = 1; %是否使用预设初始形状，0：默认均匀初始值; 1：使用给定初始值
V_constrain = 0; %0:等式体积约束，1：不等式体积约束
xinitial = 1;
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
end   
x(Ind) = 1;
if continuation == 1
    figure; imshow(1-x);
end
[y,loop,loop_k,c,x,energies,energies_k]=topthr_penlity(nelx,nely,volfrac,lambda,r,g,sd,bc,continuation,x,fileID,V_constrain);
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
