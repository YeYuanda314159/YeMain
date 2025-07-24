function [ce,cq,c] = solver_heat(xPhys,nelx,nely,freedofs)
%% MATERIAL PROPERTIES
kapa = [10, 1];  %设置热导率
q1 = 1; q2 = 100;
q    = [q1, q2]/nelx/nely; %设置热源，做尺度变换
%% PREPARE FINITE ELEMENT ANALYSIS
KE = 1/6*[4, -1, -2, -1;...
          -1, 4, -1, -2;...
          -2, -1, 4, -1;...
          -1, -2, -1, 4];    
ME = 1/4*ones(4,1);
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx); %按列取数据排成(1+nely)行,(1+nelx)列
edofVec = reshape(nodenrs(2:end,1:end-1), nelx*nely, 1); %按列取数据排成(nelx*nely)行,1列
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1], nelx*nely, 1); %给出每个元的节点的编号
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);%遍历每个单元依次经过的节点编号
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);%用于装配刚度矩阵
iM = reshape(edofMat',4*nelx*nely,1);

%% FE-ANALYSIS
%每个物理单元的刚度矩阵
sK = reshape(KE(:)*(kapa(2)+xPhys(:)'*(kapa(1)-kapa(2))), 16*nelx*nely, 1);
sF  = reshape(ME(:)*((q(1)-q(2))*xPhys(:)'+ q(2)), 4*nelx*nely, 1);
%重排前的矩阵每一列表示一个物理单元的刚度矩阵
%xPhys(:),表示按列取数据将矩阵变成一列向量
U = zeros((nely + 1)*(nelx + 1), 1); %解向量
K = sparse(iK,jK,sK); K = (K+K')/2; F = sparse(iM,1,sF);%装配整体刚度矩阵并保证对称性
U(freedofs) = K(freedofs,freedofs)\F(freedofs); %U是位移向量

%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
%U(edofMat)每行表示一个单元的自由节点位移
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); %每个单元的单位柔度  
cq = reshape(sum((U(edofMat)/4).*U(edofMat),2),nely,nelx);
c = sum(sum((kapa(2)+xPhys*(kapa(1)-kapa(2))).*ce)); %计算总能
end