%%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** 最小单元尺度固定为1
function [y, loop, c,  x, energies, energiessame]=topthr(nelx, nely, volfrac, frac, g, sd, bc, continuation, x, fileID, Vforce, sameEmin)
% nelx: number of elements on x axis
% nely: number of elements on y axis
% volfrac: volume fraction of material to total area
% frac: fraction of Emin to E0
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
E0 = 5000*8/3;  %设置杨氏模量
Emin = E0*frac; %设置人工材料的杨氏模量
nu = 1/3;       %设置泊松比
tol = 0.01;   %容许变化
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = (1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])); 
% Ae = 2\mu e + \lambda tr(e) I
% \mu = E /2/ (1+nu) \lambda = E*nu/(1+nu)/(1+nu*(1-d))
%matrix_dim = 8*8
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx); %按列取数据排成(1+nely)行,(1+nelx)列
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1); %按列取数据排成(nelx*mely)行,1列
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1); %给出每个元的节点的编号
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);%遍历每个单元依次经过的节点编号
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);%用于装配刚度矩阵
%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
F = zeros(2*(nely+1)*(nelx+1),1);
switch bc
    case 'mbb'
        F(2,1) = -1; %F(2,1) = -1, 表示1号节点受到向下大小为1的点力
        fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); %左边界横向、右下角固定
    case 'cantilever_rb'
        force_length = floor(nelx/8); %1/8宽的边受力
        %下边界右侧1/8宽的边受到向下大小250的面力
        %sparse(参数1，参数2，参数3，参数4，参数5);
        %参数1填入行坐标，参数2填入列坐标，参数3填入值，参数4、5填入矩阵大小
        F(2*(nely+1)*((nelx+1-force_length):(nelx+1))) = -250/nely; 
        fixeddofs = [1 : 2*(nely+1)]; %左边界固定 
    case 'cantilever_central'
        %右边中点受向下大小为1的点力
        F(2*(nely + 1) * nelx + 2*(floor(nely/2) + 1), 1) = -1;
        fixeddofs = [1 : 2*(nely+1)];  %左边界固定
end
U = zeros(2*(nely + 1)*(nelx + 1), 1); %解向量
alldofs = [1:2*(nely+1)*(nelx+1)]; %所以自由度编号
freedofs = setdiff(alldofs,fixeddofs);  %Dirichlet 边界外的自由度编号
%% INITIALIZE ITERATION
if continuation == 0 %采用均匀密度作为初值
    x = repmat(volfrac,nely,nelx);
    xPhys = x; %matrix_dim = nely*nelx
else
    if sd > 0 %磨光因子大于零，可以进行磨光
        xPhys = imgaussfilt(x, sd); %基于内部启发式方法，在空间域或频域中执行卷积
        %重复边界元素填充矩阵
    else
        xPhys = x;
    end
end
M = floor(nelx*nely*volfrac); %向下取整，\Omega1的元的数量
loop = 0;
change = 100; %两次迭代中分布场xPhys的变化
gamma = g*sqrt(2*pi)/(sd*1/nely); %sd*1 = sd, 若想对分量求和应使用norm(sd,1) sd/nely = tau
energies = []; %存储每次迭代的柔度能
energiessame = []; %存储给定的sameEmin对应的柔度能
%% START ITERATION
figure('Renderer', 'painters', 'Position', [90 90 1000 nely/nelx*1000]); 
%在显示器的(90,90)位置开辟一个100(长)*(nely/nelx*100)(高)的图像窗口
print_to_file= true;
try 
    fprintf(fileID, 'Displaying\n'); %若此行执行失败则执行catch后的内容
catch err
    disp(err.message)
    disp('Now display the output in the command window.')
    print_to_file=false;
end

while change>tol %只要变化大于0.01, 迭代不停
  loop = loop + 1;%迭代次数+1
  %% FE-ANALYSIS
  %每个物理单元的刚度矩阵
  %增加体积力：YYD
  if Vforce ~= 0
      x_vector = reshape(xPhys, nelx*nely, 1);
      for i  = 1:length(x_vector)
         F(edofMat(i,2)) = F(edofMat(i,2)) - Vforce*x_vector(i)*0.25;
         F(edofMat(i,8)) = F(edofMat(i,8)) - Vforce*x_vector(i)*0.25;
         F(edofMat(i,4)) = F(edofMat(i,4)) - Vforce*x_vector(i)*0.25;
         F(edofMat(i,6)) = F(edofMat(i,6)) - Vforce*x_vector(i)*0.25;
      end 
  end
  sK = reshape(KE(:)*(1./(1/Emin+xPhys(:)'*(1/E0-1/Emin))), 64*nelx*nely, 1);
  %重排前的矩阵每一列表示一个物理单元的刚度矩阵
  %xPhys(:),表示按列取数据将矩阵变成一列向量
  K = sparse(iK,jK,sK); K = (K+K')/2; %装配整体刚度矩阵并保证对称性
  U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U是位移向量
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  %U(edofMat)每行表示一个单元的自由节点位移
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); %每个单元的单位柔度    
  c = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ce)); %计算总柔度
  energies(loop) = c; %记录总体柔度
  %% threshold dynamics--calculate c \sigma:\sigma
  csigmasigmae = 1./(1/Emin+xPhys*(1/E0-1/Emin)).^2.*ce; %柔度矩阵
  %csigmasigmae = E0*ce;
  %为什么有平方?
  %% threshold dynamics--filtering
  %\phi = (1/E0-1/Emin)G_tau*(C0\sigma:\sigma) -
  %\gamma\sqrt(pi/\tau1)G_tau*(1-2x)
  if sd > 0
    csigmasigmae = imgaussfilt(csigmasigmae, sd, 'Padding', 'symmetric');
  end
  perterm = 1-2*x;
  if sd > 0 
    perterm = imgaussfilt(perterm, sd);
  end
  %% threshold dynamics--thresholding
  phi = (1/E0 - 1/Emin)*csigmasigmae + gamma*perterm; %计算线性化函数phi
  [sortedphi,I] = sort(phi(:)); %由小到大快速排序
  thr = sortedphi(M);
  % thresholding
  if thr > 0   
    xnew = 1 - double(phi>0);%phi(i) > 0, double(phi>0)(i) ==1
  else
    xnew = zeros(nelx*nely, 1);
    xnew(I(1:M)) = 1; %最小的M个元是新的最优区域
    xnew = reshape(xnew, nely, nelx);
  end
  if nargin > 11 
     %如果输入的参数个数大于10
     %此时使用输入的sameEmin作为人造材料的杨氏模量
    sK = reshape(KE(:)*((sameEmin+xnew(:)'*(E0 - sameEmin))), 64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2; 
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); 
    energiessame(loop) = sum(sum((sameEmin+xnew*(E0-sameEmin)).*ce));   
  end
  if sd > 0
    xPhys = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	用自身的镜面反射填充图像。
  else
    xPhys = xnew;
  end
  change = max(abs(xnew(:)-x(:)));%计算更新前后的区域的无穷范数
  x = xnew;
  %% PRINT RESULTS
  if print_to_file
    fprintf(fileID, ' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change); %mean: 计算xPhys的品均值
  else
    fprintf(' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change);
  end
  %% PLOT DENSITIES 输出优化的形状
  colormap(gray); imshow(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
sK = reshape(KE(:)*(1./(1/Emin+xPhys(:)'*(1/E0-1/Emin))),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2; 
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);     
c = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ce));
loop = loop + 1;
energies(loop) = c; %计算最后的柔度
if nargin > 11 %如果输入参数个数大于10，用输入的sameEmin计算最后结果
    sK = reshape(KE(:)*((sameEmin+xnew(:)'*(E0 - sameEmin))), 64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2; 
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); 
    energiessame(loop) = sum(sum((sameEmin+xnew*(E0-sameEmin)).*ce));   
end
sK = reshape(KE(:)*(x(:)'*(E0-Emin) + Emin), 64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2; 
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%% FINAL OBJECTIVE FUNCTION WITHOUT SMOOTHING
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); 
y = sum(sum(((E0-Emin)*x+Emin).*ce));
if print_to_file
    fprintf(fileID, ' sharp interface energy: %11.6f\n', y);
else
    fprintf(' sharp interface energy: %11.8f\n', y);
end
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