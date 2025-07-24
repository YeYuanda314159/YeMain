 %%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** ��С��Ԫ�߶ȹ̶�Ϊ1
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
kapa = [10, 1];  %�����ȵ���
q1 = 1; q2 = 100;
q    = [q1, q2]/nelx/nely; %������Դ�����߶ȱ任
%% Iterative Patameters
loop = 0;
change = 100; %���ε����зֲ���xPhys�ı仯
tol = 0.01;
max_iter = 1000;
%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
switch bc
    case 'left_Dirichlet'
        dirichlet_length = floor(nely/5);
        fixeddofs = (2*dirichlet_length+1) : (3*dirichlet_length); %��߽��м�1/5�ĳ��Ⱥ���
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
alldofs = 1:(nely+1)*(nelx+1); %�������ɶȱ��
freedofs = setdiff(alldofs,fixeddofs);  %Dirichlet �߽�������ɶȱ��
%% INITIALIZE ITERATION
if continuation == 0 %���þ����ܶ���Ϊ��ֵ
    M = floor(nelx*nely*volfrac); %����ȡ����\Omega1��Ԫ������
    x = repmat(volfrac,nely,nelx);
    xPhys = x; %matrix_dim = nely*nelx
else
    M = nnz(x);%����x�ķ���Ԫ����
    if sd > 0 %ĥ�����Ӵ����㣬���Խ���ĥ��
        xPhys = imgaussfilt(x, sd); %�����ڲ�����ʽ�������ڿռ����Ƶ����ִ�о��
        %�ظ��߽�Ԫ��������
    else
        xPhys = x;
    end
end
gamma = g*sqrt(2*pi/sd*1/nely); %sd*1 = sd, ����Է������Ӧʹ��norm(sd,1) sd/nely = tau
energies = []; %�洢ÿ�ε�������
%% START ITERATION
%figure('Renderer', 'painters', 'Position', [90 90 200 nely/nelx*200]); 
%����ʾ����(90,90)λ�ÿ���һ��100(��)*(nely/nelx*100)(��)��ͼ�񴰿�
print_to_file= true;
fid = fopen(fileID, 'wt');  % 'a' ��ʾ append ģʽ
try 
    fprintf(fid, 'Displaying\n'); %������ִ��ʧ����ִ��catch�������
    fprintf(fid,'Kapa:[%5.3f, %5.3f]|q:[%5.3f, %5.3f]|alpha:%5.6f|gamma:%5.6f\n',...
        kapa(1),kapa(2),q1,q2,lambda,gamma);
    fprintf(fid,'mesh:%5.6f | boundary condition:%s | continuation:%d\n',...
        1/nelx,bc,continuation);
catch err
    disp(err.message)   
    disp('Now display the output in the command window.')
    print_to_file=false;
end
while change>tol && loop <= max_iter %ֻҪ�仯����0.01, ������ͣ
  loop = loop + 1;%��������+1
  [ce,cq,c] = solver_heat(xPhys,nelx,nely,freedofs);
  energies(loop) = c; %��¼����

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
  [sortedphi,I] = sort(phi(:)); %��С�����������
  thr = sortedphi(M);
  %thresholding
  if thr > 0   
    xnew = 1 - double(phi>0);%phi(i) > 0, double(phi>0)(i) ==1
  else
    xnew = zeros(nelx*nely, 1);
    xnew(I(1:M)) = 1; %��С��M��Ԫ���µ���������
    xnew = reshape(xnew, nely, nelx);
  end
  if sd > 0
    xPhys = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	������ľ��淴�����ͼ��
  else
    xPhys = xnew;
  end
  change = norm(xnew-x,1);%�������ǰ�������������
  x = xnew;
  %% PRINT RESULTS
  if print_to_file
    fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f\n', loop, c, ...
    mean(xPhys(:)),change); %mean: ����xPhys��Ʒ��ֵ
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
energies(loop) = c; %����������

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