function [ce,cq,c] = solver_heat(xPhys,nelx,nely,freedofs)
%% MATERIAL PROPERTIES
kapa = [10, 1];  %�����ȵ���
q1 = 1; q2 = 100;
q    = [q1, q2]/nelx/nely; %������Դ�����߶ȱ任
%% PREPARE FINITE ELEMENT ANALYSIS
KE = 1/6*[4, -1, -2, -1;...
          -1, 4, -1, -2;...
          -2, -1, 4, -1;...
          -1, -2, -1, 4];    
ME = 1/4*ones(4,1);
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx); %����ȡ�����ų�(1+nely)��,(1+nelx)��
edofVec = reshape(nodenrs(2:end,1:end-1), nelx*nely, 1); %����ȡ�����ų�(nelx*nely)��,1��
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1], nelx*nely, 1); %����ÿ��Ԫ�Ľڵ�ı��
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);%����ÿ����Ԫ���ξ����Ľڵ���
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);%����װ��նȾ���
iM = reshape(edofMat',4*nelx*nely,1);

%% FE-ANALYSIS
%ÿ������Ԫ�ĸնȾ���
sK = reshape(KE(:)*(kapa(2)+xPhys(:)'*(kapa(1)-kapa(2))), 16*nelx*nely, 1);
sF  = reshape(ME(:)*((q(1)-q(2))*xPhys(:)'+ q(2)), 4*nelx*nely, 1);
%����ǰ�ľ���ÿһ�б�ʾһ������Ԫ�ĸնȾ���
%xPhys(:),��ʾ����ȡ���ݽ�������һ������
U = zeros((nely + 1)*(nelx + 1), 1); %������
K = sparse(iK,jK,sK); K = (K+K')/2; F = sparse(iM,1,sF);%װ������նȾ��󲢱�֤�Գ���
U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U��λ������

%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
%U(edofMat)ÿ�б�ʾһ����Ԫ�����ɽڵ�λ��
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); %ÿ����Ԫ�ĵ�λ���  
cq = reshape(sum((U(edofMat)/4).*U(edofMat),2),nely,nelx);
c = sum(sum((kapa(2)+xPhys*(kapa(1)-kapa(2))).*ce)); %��������
end