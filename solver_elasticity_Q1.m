function [UV,c] = solver_elasticity_Q1(xPhys,F,H,freedofs,nelx,nely,E0,Emin)
%% ��λ��Ԫ�նȾ���
nu = 1/3;       %���ò��ɱ�
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = (1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])); 
%% PREPARE FINITE ELEMENT ANALYSIS
% Ae = 2\mu e + \lambda tr(e) I
% \mu = E /2/ (1+nu) \lambda = E*nu/(1+nu)/(1+nu*(1-d))
%matrix_dim = 8*8
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx); %����ȡ�����ų�(1+nely)��,(1+nelx)��
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1); %����ȡ�����ų�(nelx*mely)��,1��
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1); %����ÿ��Ԫ�Ľڵ�ı��
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);%����ÿ����Ԫ���ξ����Ľڵ���
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);%����װ��նȾ���
%% ���
U = zeros(2*(nely + 1)*(nelx + 1), 1); %Լ���Ľ�����
V = zeros(2*(nely + 1)*(nelx + 1), 1); %Ŀ�꺯����Ӧ�Ľ�����
sK = reshape(KE(:)*(1./(1/Emin+xPhys(:)'*(1/E0-1/Emin))), 64*nelx*nely, 1);%�������
%����ǰ�ľ���ÿһ�б�ʾһ������Ԫ�ĸնȾ���
%xPhys(:),��ʾ����ȡ���ݽ�������һ������
K = sparse(iK,jK,sK); K = (K+K')/2; %װ������նȾ��󲢱�֤�Գ���
U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U��״̬��λ������
V(freedofs) = K(freedofs,freedofs)\H(freedofs); %V�ǰ����λ������

%% ������غ����ͷ���ֵ
%U(edofMat)ÿ�б�ʾһ����Ԫ�����ɽڵ�λ��
%ceU = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); %ÿ����Ԫ�ĵ�λ��� 
%ceV = reshape(sum((V(edofMat)*KE).*V(edofMat),2),nely,nelx);
ceW = reshape(sum((U(edofMat)*KE).*V(edofMat),2),nely,nelx);
%UU = 1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceU; %��Ⱥ���: h(\chi)C\sigma:\sigma
UV = 1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceW; %��Ⱥ���: h(\chi)C\tau:\tau
%VV = 1./(1/Emin+xPhys*(1/E0-1/Emin)).*ceV;%��Ⱥ���: h(\chi)C\sigma:\tau
%cU = sum(sum(UU)); %���������
%cV = sum(sum(VV));
c = sum(sum(UV));

end