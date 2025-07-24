%%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** ��С��Ԫ�߶ȹ̶�Ϊ1
function [y, loop, c,  x, energies, energies_k]=topthr_penlity(nelx, nely, lambda, volfrac, frac, g, sd, objectfunc, bc, w,continuation, x, fileID,logtype, Vforce)
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
E0 = 5000*8/3;  %��������ģ��
Emin = E0*frac; %�����˹����ϵ�����ģ��
nu = 1/3;       %���ò��ɱ�
tol = 0.01;   %����仯
p   =  1;     %��ʽ������
maxtimes = 100;  %����������
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = (1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])); 
% Ae = 2\mu e + \lambda tr(e) I
% \mu = E /2/ (1+nu) \lambda = E*nu/(1+nu)/(1+nu*(1-d))
%matrix_dim = 8*8
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx); %����ȡ�����ų�(1+nely)��,(1+nelx)��
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1); %����ȡ�����ų�(nelx*mely)��,1��
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1); %����ÿ��Ԫ�Ľڵ�ı��
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);%����ÿ����Ԫ���ξ����Ľڵ���
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);%����װ��նȾ���
%% DEFINE LOADS AND SUPPORTS assuming the mesh size is 1/nely
Fin = zeros(2*(nely+1)*(nelx+1),1);
Hout = zeros(2*(nely+1)*(nelx+1),1);
switch bc
    case 'left_bdc_right_up_qin'
        Fin(2*(nely+1)*nelx+2,1) = -1; %F(2,1) = -1, ��ʾ���Ͻڵ��ܵ����´�СΪ1�ĵ���
        fixeddofs = 1:2*(nely+1); %��߽����
    case 'left_down_bdc_right_down_qin'
        force_length = floor(nelx/10); %1/8��ı�����
        %sparse(����1������2������3������4������5);
        %����1���������꣬����2���������꣬����3����ֵ������4��5��������С
        Fin((2*(nelx+1)*(nely+1)-2*force_length-1):2:(2*(nelx+1)*(nely+1)-1)) = -2; 
        fixeddofs = union(1 :2: 2*(nely+1)-1,(2*(nely+1)*(1:(nelx+1)))); 
        %fixeddofs = union(fixeddofs,2*(nely+1)*(1:(nelx+1))-1);
    case 'left_up_down_bdc_left_down_qin'
        force_length = floor(nely/10);
        Fin(2*(nely + 1-force_length)+1:2:(2*(nely+1)-1)) = 1;
        fixeddofs = union(1 : 2*force_length,(2*(nely+1)*(1:(nelx+1))));  
end
switch objectfunc %���������Dirichlet�߽�ͬԪ������ͬ
    case 'down_central'
        %�±��е�ĵ�λλ��
        Hout(2*(nely + 1) * floor((nelx+1)/2), 1) = 1;
    case 'left_up'
        force_length = floor(nelx/10);
        Hout(2+(2*(nely+1)*(0:force_length))) = -1; %H(2,1) = -1, ��ʾ1�Žڵ��ܵ����´�СΪ1�ĵ���
    case 'right_down'
        force_length = floor(nely/10); %1/8��ı�����
        %�±߽��Ҳ�1/8��ı��ܵ����´�С250������
        %sparse(����1������2������3������4������5);
        %����1���������꣬����2���������꣬����3����ֵ������4��5��������С
        Hout((2*(nelx+1)*(nely+1)-2*force_length-1):2:(2*(nelx+1)*(nely+1)-1)) = 1; 
end
F = Fin + Hout;
H = Fin + w*Hout;
U = zeros(2*(nely + 1)*(nelx + 1), 1); %Լ���Ľ�����
alldofs = [1:2*(nely+1)*(nelx+1)]; %�������ɶȱ��
freedofs = setdiff(alldofs,fixeddofs);  %Dirichlet �߽�������ɶȱ��
V = zeros(2*(nely + 1)*(nelx + 1), 1); %Ŀ�꺯����Ӧ�Ľ�����
%% INITIALIZE ITERATION
if continuation == 0 %���þ����ܶ���Ϊ��ֵ
    x = repmat(volfrac,nely,nelx);
    xPhys = x.^p; %matrix_dim = nely*nelx
else
    if sd > 0 %ĥ�����Ӵ����㣬���Խ���ĥ��
        xPhys = imgaussfilt(x, sd); %�����ڲ�����ʽ�������ڿռ����Ƶ����ִ�о��
        xPhys = xPhys.^p;
        %�ظ��߽�Ԫ��������
    else
        xPhys = x;
    end
end
M = floor(nelx*nely*volfrac); %����ȡ����\Omega1��Ԫ������
loop = 0;
change = 1; %���ε����зֲ���xPhys�ı仯
gamma = g*sqrt(2*pi)/(sd*1/nely); %sd*1 = sd, ����Է������Ӧʹ��norm(sd,1) sd/nely = tau
energies = []; %�洢ÿ�ε����������
energies_k = [];
%% START ITERATION
figure('Renderer', 'painters', 'Position', [90 90 1000 nely/nelx*1000]); 
%����ʾ����(90,90)λ�ÿ���һ��100(��)*(nely/nelx*100)(��)��ͼ�񴰿�
fid = fopen(fileID,logtype);
print_to_file= true;
try 
    fprintf(fid, 'Displaying\n'); %������ִ��ʧ����ִ��catch�������
    fprintf(fid, 'bc:%s | objectfunc:%s | Vforce:%3.5f | delta_x:%3.5f | volfrac:%3.5f\n',...
        bc,objectfunc,Vforce,1/nely,volfrac);
    fprintf(fid, 'Emin:%5.5f | E0:%5.5f | gamma:%3f | tau:%5.5f | lambda:%5.5f | w:%5.5f\n',...
        Emin, E0, gamma, sd/nely, lambda,w);
catch err
    disp(err.message)
    disp('Now display the output in the command window.')
    print_to_file=false;
end

while change>tol && loop < maxtimes%ֻҪ�仯����0.01, ������ͣ, ������maxtimes��
  loop = loop + 1;%��������+1
  %% FE-ANALYSIS
  %ÿ������Ԫ�ĸնȾ���
  %�����������YYD
  if Vforce ~= 0
      x_vector = reshape(xPhys, nelx*nely, 1);
      for i  = 1:length(x_vector)
         F(edofMat(i,2)) = F(edofMat(i,2)) - Vforce*x_vector(i)*0.25;
         F(edofMat(i,8)) = F(edofMat(i,8)) - Vforce*x_vector(i)*0.25;
         F(edofMat(i,4)) = F(edofMat(i,4)) - Vforce*x_vector(i)*0.25;
         F(edofMat(i,6)) = F(edofMat(i,6)) - Vforce*x_vector(i)*0.25;
      end 
  end
  sK = reshape(KE(:)*(1./(1/Emin+xPhys(:)'*(1/E0-1/Emin))), 64*nelx*nely, 1);%�������
  %����ǰ�ľ���ÿһ�б�ʾһ������Ԫ�ĸնȾ���
  %xPhys(:),��ʾ����ȡ���ݽ�������һ������
  K = sparse(iK,jK,sK); K = (K+K')/2; %װ������նȾ��󲢱�֤�Գ���
  U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U��״̬��λ������
  V(freedofs) = K(freedofs,freedofs)\H(freedofs); %V�ǰ����λ������
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  %U(edofMat)ÿ�б�ʾһ����Ԫ�����ɽڵ�λ��
  ceU = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); %ÿ����Ԫ�ĵ�λ��� 
  ceV = reshape(sum((V(edofMat)*KE).*V(edofMat),2),nely,nelx);
  ceW = reshape(sum((U(edofMat)*KE).*V(edofMat),2),nely,nelx);
  cU = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceU)); %���������
  cV = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceV));
  c = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceW));
  energies(loop) = c; %��¼�������
  %% threshold dynamics--calculate c \sigma:\sigma
  perterm = 1-2*x;
  if sd > 0 
    perterm = imgaussfilt(perterm, sd);
  end
  csigmasigma = 1./(1/Emin+xPhys*(1/E0-1/Emin)).^2.*(ceW);%  + lambda*(ceU+ceV)); %��Ⱦ���
  
  %% threshold dynamics--filtering
  %\phi = (1/E0-1/Emin)G_tau*(C0\sigma:\sigma) -
  %\gamma\sqrt(pi/\tau1)G_tau*(1-2x)
  %% threshold dynamics--thresholding
  if sd > 0
    csigmasigma = imgaussfilt(csigmasigma, sd, 'Padding', 'symmetric');
  end
  phi = (1/E0 - 1/Emin)*csigmasigma + gamma*perterm; %�������Ի�����phi
  [sortedphi,I] = sort(phi(:)); %��С�����������
  thr = sortedphi(M);
  % thresholding
  if thr > 0   
    xnew = 1 - double(phi>0);%phi(i) > 0, double(phi>0)(i) ==1
  else
    xnew = zeros(nelx*nely, 1);
    xnew(I(1:M)) = 1; %��С��M��Ԫ���µ���������
    xnew = reshape(xnew, nely, nelx);
  end
  if sd > 0
    xPhys = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	������ľ��淴�����ͼ��
    xPhys = xPhys.^p;
  else
    xPhys = xnew.^p;
  end
  change = norm(xnew-x,1);%�������ǰ�������������
  x = xnew;
  %����������º�ĺ���ֵ
  cUPhy = sum(sum((1/Emin+xPhys*(1/E0-1/Emin)).*ceU)); %���������
  cVPhy = sum(sum((1/Emin+xPhys*(1/E0-1/Emin)).*ceV)); 
  cWPhy = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceW));
  cPhy =  cWPhy + lambda*(cUPhy + cVPhy - cU - cV);
  energies_k(loop) = cPhy; %��¼�������
  %% PRINT RESULTS
  if print_to_file
    fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||objxPhy.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f\n', loop, c, ...
    cPhy, mean(xPhys(:)),change); %mean: ����xPhys��Ʒ��ֵ
  else
    fprintf(' It.:%5i Obj.:%10.6f objxPhy.:%10.6f Vol.:%7.3f ch.:%7.3f\n', loop, c, cPhy,...
    mean(xPhys(:)),change);
  end
  %% PLOT DENSITIES ����Ż�����״
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
sK = reshape(KE(:)*(1./(1/Emin+xPhys(:)'*(1/E0-1/Emin))), 64*nelx*nely, 1);%�������
  %����ǰ�ľ���ÿһ�б�ʾһ������Ԫ�ĸնȾ���
  %xPhys(:),��ʾ����ȡ���ݽ�������һ������
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U��״̬��λ������
V(freedofs) = K(freedofs,freedofs)\H(freedofs); %V�ǰ����λ������
ceW = reshape(sum((U(edofMat)*KE).*V(edofMat),2),nely,nelx);
c = sum(sum(1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceW));
loop = loop + 1;
energies(loop) = c; %�����������
sK = reshape(KE(:)*(x(:)'*(E0-Emin) + Emin), 64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2; 
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
V(freedofs) = K(freedofs,freedofs)\H(freedofs);
%% FINAL OBJECTIVE FUNCTION WITHOUT SMOOTHING
ceW = reshape(sum((U(edofMat)*KE).*V(edofMat),2),nely,nelx);
y = sum(sum(((E0-Emin)*x+Emin).*ceW));
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