%%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** 最小单元尺度固定为1
function [y, loop, loop_k, c,  x, energies, energies_k]=topthr_penlity(nelx, nely, volfrac, lambda, r0, g, sd, bc, continuation, x, fileID, V_constrain)
% nelx: number of elements on x axis
% nely: number of elements on y axis
% volfrac: volume fraction of material to total area
% alpha: coefficient of penalty term
% sd1: filter size for thresholding
% sd2: filter size for solving equation
% sd3: filter size for the perimeter term
% bc:  takes values in {'cantilever_rb', 'cantilever_central', 'mbb'},
%      enforces BCs (boundary condition)
% g: coefficient of the perimeter term
% continuation: 1 indicates to use x as the initial guess
%               0 indicates to use constant density as the initial guess
% x: initial guess, if continuaton is 0, can be anything.
% fileID: the opened file to log outputs
% If you don't want to save the results, just let fileID to be negative 
% and the results will be displayed in the command window
% sameEmin: calculate the compliance each iteration using this specified value as Emin
% return c: compliance computed with filtered chi
%% MATERIAL PROPERTIES
kapa = [10, 1];  %设置热导率
q1 = 1; q2 = 100;
q    = [q1, q2]/nelx/nely; %设置热源，做尺度变换
%% Iterative Patameters
loop = 1; %变量迭代循环
loop_k = 1; %线搜索循环
change = 100; %两次迭代中分布场xPhys的变化
tol = 2;
max_iter = 500;
r = r0; %默认最大邻近系数
gamma = g*sqrt(2*pi/sd*1/nely); %sd*1 = sd, 若想对分量求和应使用norm(sd,1) sd/nely = tau
energies = []; %存储每次迭代的能
energies_k = [];
%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
switch bc
    case 'left_Dirichlet'
        dirichlet_length = floor(nely/5);
        fixeddofs = (2*dirichlet_length+1) : (3*dirichlet_length); %左边界中间1/5的长度恒温
    case 'all_Dirichlet'
        rig = 1:(nely+1);
        top = (nely+2):(nely+1):((nelx-1)*(nely+1)+1);
        riglef = union(rig,(nely+1)*nelx + rig);
        uppdow = union(top,top+nely);
        fixeddofs = union(riglef,uppdow);
    case 'topleft_Dirichlet'
        rig = 1:(nely+1);
        top = (nely+2):(nely+1):((nelx-1)*(nely+1)+1);
        fixeddofs = union(rig,top);
    case 'allleft_Dirichlet'
        fixeddofs = 1:(nely+1);
end
alldofs = 1:(nely+1)*(nelx+1); %所有自由度编号
freedofs = setdiff(alldofs,fixeddofs);  %Dirichlet 边界外的自由度编号
%% INITIALIZE ITERATION
if continuation == 0 %采用均匀密度作为初值
    M = floor(nelx*nely*volfrac); %向下取整，\Omega1的元的数量
    x = repmat(volfrac,nely,nelx);
    xPhys = x; %matrix_dim = nely*nelx
else
    M = nnz(x);%计算x的非零元个数
    if sd > 0 %磨光因子大于零，可以进行磨光
        xPhys = imgaussfilt(x, sd); %基于内部启发式方法，在空间域或频域中执行卷积
        %重复边界元素填充矩阵
    else
        xPhys = x;
    end
end
%% START ITERATION
%figure('Renderer', 'painters', 'Position', [90 90 200 nely/nelx*200]); 
%在显示器的(90,90)位置开辟一个100(长)*(nely/nelx*100)(高)的图像窗口
print_to_file= true;
fid = fopen(fileID, 'wt');  % 'a' 表示 append 模式
try 
    fprintf(fid, 'Displaying\n'); %若此行执行失败则执行catch后的内容
    fprintf(fid,'Kapa:[%5.3f, %5.3f]|q:[%5.3f, %5.3f]|lambda:%5.6f | r:%5.6f |gamma:%5.6f\n',...
        kapa(1),kapa(2),q1,q2,lambda,r,gamma);
    fprintf(fid,'mesh:%5.6f | boundary condition:%s | continuation:%d\n',...
        1/nelx,bc,continuation);
catch err
    disp(err.message)   
    disp('Now display the output in the command window.')
    print_to_file=false;
end
%创建固定的图像窗口
ax = axes('Parent', gcf);
while 1
    if change <= tol  %只要变化大于0.01, 迭代不停
        fprintf('收敛原因2：自变量无法更新！\n');
        break;
    elseif loop > max_iter %最多计算max_iter次
        fprintf('收敛原因3：达到最大迭代次数！\n');
        break;
    end
    if loop == 1
        [ce,cq,c] = solver_heat(xPhys,nelx,nely,freedofs);
        PG = gamma*sum(sum((1-x).*xPhys));
        energies(loop) = c + PG; %记录总能
        energies_k(loop) = energies(loop); %记录总能
    end
    %% Penalty Method--calculate g^k = （1/(2\ambda)-1）*(kapa(1)-kapa(2))*ce + (2-1/lambda)*(q(1)-q(2))*cq
    gk = (0.5/lambda-1)*(kapa(1)-kapa(2))*ce + (2-1/lambda)*(q(1)-q(2))*cq;
    if sd > 0
        gk = imgaussfilt(gk, sd, 'Padding', 'symmetric');
    end
    gk = (gk + gamma*(x-xPhys))/(1.5-1/lambda);
    if loop == 1
        gk0 = gk;
        wk = gk;
    else
        betak = max(sum(sum(gk.*(gk-gk0)))/norm(gk0,2)^2,0); 
        wk = gk + betak*wk;
        gk0 = gk;
    end
    phi = wk;
    %% Penalty Method--linear research parameter
    if loop == 1
        r_min = 0;
    else
        %找出x中取值为1的索引
        mask = x > 0;
        indx1 = find(mask);
        indx0 = setdiff(1:numel(x), indx1); % 索引0（补集）
        % 最高效方法：逻辑索引
        max_val_indx1 = max(phi(indx1));  % indx1中最大值
        min_val_indx0 = min(phi(indx0));  % indx0中最小值
        if max_val_indx1 <= min_val_indx0 %满足一阶最优性条件，迭代终止
            fprintf('收敛原因1：达到一阶最优性条件！\n');
            break;
        end
        r_min = 1/(max_val_indx1 - min_val_indx0);
    end
    sorted_A = sort(phi(:));  % 从大到小排序成向量
    min_spacing = min(diff(sorted_A)); % 相邻元素最小间距
    if r <= r_min + 1;
        r = min(1/min_spacing,r0);
    end
    %% Penalty Method--linear research
    while 1
        if (loop_k - loop) > max_iter
            change = 0;
            fprintf('线搜索终止3：达到最大线搜索次数！\n');
            break;
        end
        loop_k = loop_k + 1;
        bar_phi = x - r*phi;   
        if sd > 0
            bar_phi = imgaussfilt(bar_phi, sd, 'Padding', 'symmetric');
        end
        [~,I] = sort(bar_phi(:),'descend'); %由大到小快速排序
        % Project
        xnew = zeros(nelx*nely, 1);
        if (bar_phi(I(M)) > 0.5) || (V_constrain == 0)
            xnew(I(1:M)) = 1; %最大的M个元是新的最优区域
        else
            xnew(bar_phi > 0.5) = 1;
        end
        xnew = reshape(xnew, nely, nelx);
        if sd > 0
            xnewPhys = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	用自身的镜面反射填充图像。
        else
            xnewPhys = xnew;
        end
        %% PLOT DENSITIES
        [ce,cq,c] = solver_heat(xnewPhys,nelx,nely,freedofs);
        til_c = c + gamma*sum(sum((1-x).*xPhys));
        energies_k(loop_k) = til_c;
        change = norm(xnew-x,1);%计算更新前后的区域的无穷范数
        %输出新图像
        cla(ax);  % 清除当前axes内容（保留axes设置）
        imshow(1-xnew, [], 'Parent', ax);  % 在同一axes显示新图像
        title(ax, sprintf('迭代 %d/%d', loop_k, max_iter));
        colormap(ax, gray);
        drawnow;
        %% PRINT RESULTS
        if print_to_file
            fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f || r.%5.9f \n', loop, c, ...
            mean(xPhys(:)),change,r); %mean: 计算xPhys的品均值
        else
            fprintf(' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f\n', loop, c, ...
            mean(xPhys(:)),change);
        end
        if til_c > energies(loop)
          r_max = r;
          if change <= 1
              change = 0;
              fprintf('线搜索终止1：自变量无法更新！\n');
              break;
          end
        elseif (til_c < energies(loop)) 
          x = xnew;
          xPhys = xnewPhys;
          loop = loop + 1;
          energies(loop) = til_c;
          fprintf('线搜索终止2：目标函数下降！\n');
          break;
        elseif (til_c == energies(loop)) && (change > tol)
          r_max = r;
        end
        if change == 0;
          r_min = r;
        end
        r = (r_max + r_min)/2;
        if abs(r_max-r_min) < 10^-3
            change = 0;
            fprintf('线搜索终止1：自变量无法更新！\n');
            break;
        end
    end
end
set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
[~,~,c] = solver_heat(xPhys,nelx,nely,freedofs);
loop = loop + 1;
energies(loop) = c + gamma*sum(sum((1-x).*xPhys)); %计算最后的能

%% FINAL OBJECTIVE FUNCTION WITHOUT SMOOTHING
[~,~,y] = solver_heat(x,nelx,nely,freedofs);
if print_to_file
    fprintf(fid, 'Thermal dissipation energy: %11.6f \n', y);
    fprintf(fid, '-----------------------------------------\n');
    fprintf(fid, '-----------------------------------------\n');
else
    fprintf(' Thermal dissipation energy: %11.8f \n', y);
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written based on 88 line MATLAB code written by 
% E. Andreassen, A. Clausen, M. Schevenels,
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
% Technical University of Denmark,                                         %
% DK-2800 Lyngby, Denmark.                                                 %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%