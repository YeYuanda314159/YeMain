%%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** 最小单元尺度固定为1
function [y, loop, c,  x, energies, energies_k]=topthr_penlity_general(nelx, nely, lambda, r, volfrac, frac, g, sd, objectfunc, bc, w,continuation, x, fileID,logtype, Vforce)
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
%% OPTIMIZATION PARAMETERS
tol = 0.01;   %容许变化
maxtimes = 100;  %最大迭代次数
gamma1 = 3;
gamma2 = 1.5;
%% MATERIAL PROPERTIES
E0 = 5000*8/3;  %设置杨氏模量
Emin = E0*frac; %设置人工材料的杨氏模量

%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
Fin = zeros(2*(nely+1)*(nelx+1),1);
Hout = zeros(2*(nely+1)*(nelx+1),1);
switch bc
    case 'left_bdc_right_up_qin'
        Fin(2*(nely+1)*nelx+2,1) = -1; %F(2,1) = -1, 表示右上节点受到向下大小为1的点力
        fixeddofs = 1:2*(nely+1); %左边界横向
    case 'left_down_bdc_right_down_qin'
        force_length = floor(nelx/10); %1/8宽的边受力
        %sparse(参数1，参数2，参数3，参数4，参数5);
        %参数1填入行坐标，参数2填入列坐标，参数3填入值，参数4、5填入矩阵大小
        Fin((2*(nelx+1)*(nely+1)-2*force_length-1):2:(2*(nelx+1)*(nely+1)-1)) = -2; 
        fixeddofs = union(1 :2: 2*(nely+1)-1,(2*(nely+1)*(1:(nelx+1)))); 
        %fixeddofs = union(fixeddofs,2*(nely+1)*(1:(nelx+1))-1);
    case 'left_up_down_bdc_left_down_qin'
        force_length = floor(nely/10);
        Fin(2*(nely + 1-force_length)+1:2:(2*(nely+1)-1)) = 1;
        fixeddofs = union(1 : 2*force_length,(2*(nely+1)*(1:(nelx+1))));  
end
switch objectfunc %伴随问题的Dirichlet边界同元问题相同
    case 'down_central'
        %下边中点的单位位移
        Hout(2*(nely + 1) * floor((nelx+1)/2), 1) = 1;
    case 'left_up'
        force_length = floor(nelx/10);
        Hout(2+(2*(nely+1)*(0:force_length))) = -1; %H(2,1) = -1, 表示1号节点受到向下大小为1的点力
    case 'right_down'
        force_length = floor(nely/10); %1/8宽的边受力
        %下边界右侧1/8宽的边受到向下大小250的面力
        %sparse(参数1，参数2，参数3，参数4，参数5);
        %参数1填入行坐标，参数2填入列坐标，参数3填入值，参数4、5填入矩阵大小
        Hout((2*(nelx+1)*(nely+1)-2*force_length-1):2:(2*(nelx+1)*(nely+1)-1)) = 1; 
end
F = Fin + Hout;
H = Fin + w*Hout;
alldofs = [1:2*(nely+1)*(nelx+1)]; %所有自由度编号
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
change = 1; %两次迭代中分布场xPhys的变化
gamma = g*sqrt(2*pi)/(sd*1/nely); %sd*1 = sd, 若想对分量求和应使用norm(sd,1) sd/nely = tau
energies = []; %存储每次迭代的柔度能
energies_k = [];
%% START ITERATION
%figure('Renderer', 'painters', 'Position', [90 90 1000 nely/nelx*1000]); 
%在显示器的(90,90)位置开辟一个100(长)*(nely/nelx*100)(高)的图像窗口
fid = fopen(fileID,logtype);
print_to_file= true;
try 
    fprintf(fid, 'Displaying\n'); %若此行执行失败则执行catch后的内容
    fprintf(fid, 'bc:%s | objectfunc:%s | Vforce:%3.5f | delta_x:%3.5f | volfrac:%3.5f\n',...
        bc,objectfunc,Vforce,1/nely,volfrac);
    fprintf(fid, 'Emin:%5.5f | E0:%5.5f | gamma:%3f | tau:%5.5f | lambda:%5.5f | w:%5.5f | gamma1:%5.5f | gamma2:%5.5f\n',...
        Emin, E0, gamma, sd/nely, lambda,w, gamma1,gamma2);
catch err
    disp(err.message)
    disp('Now display the output in the command window.')
    print_to_file=false;
end

while change>tol && loop < maxtimes%只要变化大于0.01, 迭代不停, 最多计算maxtimes次
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
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [UV,c] = solver_elasticity_Q1(xPhys,F,H,freedofs,nelx,nely,E0,Emin);
  energies(loop) = c; %记录总体柔度
 
  %% threshold dynamics--calculate H = diag{aa+bb} - 1/|aa| * aa*aa' - 1/|bb| * bb*bb'
  % D = diag{aa+bb}, B = [1/|a|*(aa),1/|b|*(bb)] => H = D - B*B'
  % H^-1 = D^-1 + D^-1*B(I-B'*D^-1*B)^-1*B'*D^-1
  % Phi = H^-1*ab = (D^-1 + D^-1*B*(I-B'*D^-1*B)^-1*B'*D^-1) * ab
  %               = ab./(aa+bb) + B./(aa+bb)*(I-B'*B./(aa+bb))^-1*B'*ab./(aa+bb)
%   aa = reshape(UU, nelx*nely, 1);
%   bb = reshape(VV, nelx*nely, 1);
%   ab = reshape(UV, nelx*nely, 1);
%   a = sqrt(sum(aa)); b = sqrt(sum(bb));
%   cc = aa + bb + (xPhys(:)-E0/(E0-Emin)).^2/r;
%   ab_bar = ab./cc;
%   B = [aa/a,bb/b]; B_bar = [B(:,1)./cc, B(:,2)./cc];
%   Phi = ab_bar + B_bar*((eye(2)-B'*B_bar)\B')*ab_bar;
%   Phi = (xPhys(:)-E0/(E0-Emin))/2/lambda.*Phi;
  Phi = 1./(1/Emin+xPhys*(1/E0-1/Emin)).*UV;
  if sd > 0 
    Phi = imgaussfilt(Phi, sd, 'Padding', 'symmetric');
  end 
 %% threshold dynamics--filtering
  perterm = 1-2*x;
  
  %\gamma\sqrt(pi/\tau1)G_tau*(1-2x)
  %% threshold dynamics--thresholding
  Phi = x - r*(1/E0-1/Emin)*Phi; %计算修正局部极小Phi
  if sd > 0 
    perterm = imgaussfilt(perterm, sd);
    Phi = imgaussfilt(Phi, sd, 'Padding', 'symmetric') + g*perterm;
  end  
  [~,I] = sort(Phi(:),'descend'); %由大到小快速排序
  
  % thresholding
  xnew = zeros(nelx*nely, 1);
  xnew(I(1:M)) = 1; %最大的M个元是新的最优区域
  xnew = reshape(xnew, nely, nelx);
  %调整邻近因子
  if sd > 0
    xPhys_new = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	用自身的镜面反射填充图像。
  else
    xPhys_new = xnew;
  end
  [~,c] = solver_elasticity_Q1(xPhys_new,F,H,freedofs,nelx,nely,E0,Emin);
  if (energies(loop)-c)/energies(loop) < 0.01 &&  energies(loop)>c
      r = r*gamma1;
  elseif energies(loop)<c
      r = r/gamma2;
      continue;
  end
  change = norm(xnew-x,1);%计算更新前后的区域的无穷范数
  x = xnew;
  xPhys = xPhys_new;  
  %% PRINT RESULTS
  if print_to_file
    fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f ||r.:%7.3f \n', loop, c, ...
    mean(xPhys(:)),change,r); %mean: 计算xPhys的品均值
  else
    fprintf(' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f r.%7.3f:\n', loop, c,...
    mean(xPhys(:)),change,r);
  end
  %% PLOT DENSITIES 输出优化的形状
  colormap(gray); imshow(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
set(gca,'Units','normalized','Position',[0 0 10 10]);  %# Modify axes size
if Vforce ~= 0
  x_vector = reshape(xPhys, nelx*nely, 1);
  for i  = 1:length(x_vector)
     F(edofMat(i,2)) = F(edofMat(i,2)) - Vforce*x_vector(i)*0.25;
     F(edofMat(i,8)) = F(edofMat(i,8)) - Vforce*x_vector(i)*0.25;
     F(edofMat(i,4)) = F(edofMat(i,4)) - Vforce*x_vector(i)*0.25;
     F(edofMat(i,6)) = F(edofMat(i,6)) - Vforce*x_vector(i)*0.25;
  end 
end
[~,c] = solver_elasticity_Q1(xPhys,F,H,freedofs,nelx,nely,E0,Emin);
loop = loop + 1;
energies(loop) = c; %计算最后的柔度
%% FINAL OBJECTIVE FUNCTION WITHOUT SMOOTHING
[~,y] = solver_elasticity_Q1(x,F,H,freedofs,nelx,nely,E0,Emin);
if print_to_file
    fprintf(fid, ' sharp interface energy: %11.6f\n', y);
    fprintf(fid, '--------------------------------------\n');
    fprintf(fid, '--------------------------------------\n');
else
    fprintf(' sharp interface energy: %11.8f\n', y);
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