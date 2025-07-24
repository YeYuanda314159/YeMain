%%%% THRESHOLD DYNAMICS USING IMGAUSSFILT %%%%
%** ��С��Ԫ�߶ȹ̶�Ϊ1
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
tol = 0.01;   %����仯
maxtimes = 100;  %����������
gamma1 = 3;
gamma2 = 1.5;
%% MATERIAL PROPERTIES
E0 = 5000*8/3;  %��������ģ��
Emin = E0*frac; %�����˹����ϵ�����ģ��

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
alldofs = [1:2*(nely+1)*(nelx+1)]; %�������ɶȱ��
freedofs = setdiff(alldofs,fixeddofs);  %Dirichlet �߽�������ɶȱ��
%% INITIALIZE ITERATION
if continuation == 0 %���þ����ܶ���Ϊ��ֵ
    x = repmat(volfrac,nely,nelx);
    xPhys = x; %matrix_dim = nely*nelx
else
    if sd > 0 %ĥ�����Ӵ����㣬���Խ���ĥ��
        xPhys = imgaussfilt(x, sd); %�����ڲ�����ʽ�������ڿռ����Ƶ����ִ�о��
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
%figure('Renderer', 'painters', 'Position', [90 90 1000 nely/nelx*1000]); 
%����ʾ����(90,90)λ�ÿ���һ��100(��)*(nely/nelx*100)(��)��ͼ�񴰿�
fid = fopen(fileID,logtype);
print_to_file= true;
try 
    fprintf(fid, 'Displaying\n'); %������ִ��ʧ����ִ��catch�������
    fprintf(fid, 'bc:%s | objectfunc:%s | Vforce:%3.5f | delta_x:%3.5f | volfrac:%3.5f\n',...
        bc,objectfunc,Vforce,1/nely,volfrac);
    fprintf(fid, 'Emin:%5.5f | E0:%5.5f | gamma:%3f | tau:%5.5f | lambda:%5.5f | w:%5.5f | gamma1:%5.5f | gamma2:%5.5f\n',...
        Emin, E0, gamma, sd/nely, lambda,w, gamma1,gamma2);
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
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [UV,c] = solver_elasticity_Q1(xPhys,F,H,freedofs,nelx,nely,E0,Emin);
  energies(loop) = c; %��¼�������
 
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
  Phi = x - r*(1/E0-1/Emin)*Phi; %���������ֲ���СPhi
  if sd > 0 
    perterm = imgaussfilt(perterm, sd);
    Phi = imgaussfilt(Phi, sd, 'Padding', 'symmetric') + g*perterm;
  end  
  [~,I] = sort(Phi(:),'descend'); %�ɴ�С��������
  
  % thresholding
  xnew = zeros(nelx*nely, 1);
  xnew(I(1:M)) = 1; %����M��Ԫ���µ���������
  xnew = reshape(xnew, nely, nelx);
  %�����ڽ�����
  if sd > 0
    xPhys_new = imgaussfilt(xnew, sd, 'Padding', 'symmetric');%	������ľ��淴�����ͼ��
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
  change = norm(xnew-x,1);%�������ǰ�������������
  x = xnew;
  xPhys = xPhys_new;  
  %% PRINT RESULTS
  if print_to_file
    fprintf(fid, ' It.:%5i ||Obj.:%10.6f ||Vol.:%7.3f ||ch.:%7.3f ||r.:%7.3f \n', loop, c, ...
    mean(xPhys(:)),change,r); %mean: ����xPhys��Ʒ��ֵ
  else
    fprintf(' It.:%5i Obj.:%10.6f Vol.:%7.3f ch.:%7.3f r.%7.3f:\n', loop, c,...
    mean(xPhys(:)),change,r);
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
[~,c] = solver_elasticity_Q1(xPhys,F,H,freedofs,nelx,nely,E0,Emin);
loop = loop + 1;
energies(loop) = c; %�����������
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