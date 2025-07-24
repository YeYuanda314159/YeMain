 %%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** 最小单元尺度固定为1
function [y, loop, c,  x, energies]=topthr_penlity(nelx, nely, volfrac, lambda, g, sd, bc, continuation, x, fileID)
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
loop = 0;
change = 100; %两次迭代中分布场xPhys的变化
tol = 0.01;
max_iter = 1000;
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
gamma = g*sqrt(2*pi/sd*1/nely); %sd*1 = sd, 若想对分量求和应使用norm(sd,1) sd/nely = tau
energies = []; %存储每次迭代的能
%% START ITERATION
%figure('Renderer', 'painters', 'Position', [90 90 200 nely/nelx*200]); 
%在显示器的(90,90)位置开辟一个100(长)*(nely/nelx*100)(高)的图像窗口
print_to_file= true;
fid = fopen(fileID, 'wt');  % 'a' 表示 append 模式
try 
    fprintf(fid, 'Displaying\n'); %若此行执行失败则执行catch后的内容
    fprintf(fid,'Kapa:[%5.3f, %5.3f]|q:[%5.3f, %5.3f]|alpha:%5.6f|gamma:%5.6f\n',...
        kapa(1),kapa(2),q1,q2,lambda,gamma);
    fprintf(fid,'mesh:%5.6f | boundary condition:%s | continuation:%d\n',...
        1/nelx,bc,continuation);
catch err
    disp(err.message)   
    disp('Now display the output in the command window.')
    print_to_file=false;
end
while change>tol && loop <= max_iter %只要变化大于0.01, 迭代不停
  loop = loop + 1;%迭代次数+1
  [ce,cq,c] = solver_heat(xPhys,nelx,nely,freedofs);
  energies(loop) = c; %记录总能

  %% threshold dynamics--filtering
  %\phi = (1/kapa(1)-1/kapa(2))G_tau*(\sigma:\sigma) - 
  %\gamma\sqrt(pi/\tau1)G_tau*(1-2x)
  if sd > 0
    ce = imgaussfilt(ce, sd, 'Padding', 'symmetric');
    cq = imgaussfilt(cq, sd, 'Padding', 'symmetric');
  end
  perterm = 1-2*x;
  if sd > 0 
    perterm = imgaussfilt(perterm, sd);
  end
  %% threshold dynamics--thresholding
  phi = -(kapa(1)-kapa(2))*ce + 2*(q(1)-q(2))*cq + gamma*perterm;
  [sortedphi,I] = sort(phi(:)); %由小到大快速排序
  thr = sortedphi(M);
  %thresholding
  if thr > 0   
    xnew = 1 - double(phi>0);%phi(i) > 0, double(phi>0)(i) ==1
  else
    xnew = zeros(nelx*nely, 1);
    xnew(I(1:M)) = 1; %最小的M个元是新的最优区域
    xnew = reshape(xnew, nely, nelx);
  end
  if sd > 0
    xPhys = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	用自身的镜面反射填充图像。
  else
    xPhys = xnew;
  end
  change = norm(xnew-x,1);%计算更新前后的区域的无穷范数
  x = xnew;
  %% PRINT RESULTS
  if print_to_file
    fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change); %mean: 计算xPhys的品均值
  else
    fprintf(' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change);
  end
  %% PLOT DENSITIES
  colormap(gray); imshow(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
[~,~,c] = solver_heat(xPhys,nelx,nely,freedofs);
loop = loop + 1;
energies(loop) = c; %计算最后的能

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