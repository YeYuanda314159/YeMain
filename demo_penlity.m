% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID = 'D:\Study\thresholding dynamics\ICTM-heat-penlity\example1.log';
% whether calculate compliance with respect to the same Emin
% parameters
nelx = 200; 
nely = 200;
volfrac = 0.2; %体积占比
lambda = 10000; %正则参数
bc = 'left_Dirichlet'; %左端1/5Dirichlet边界
%bc = 'all_Dirichlet'; %完全Dirichlet边界
%bc = 'topleft_Dirichlet';
%bc = 'allleft_Dirichlet';
continuation = 1; %是否使用预设初始形状，0：默认均匀初始值; 1：使用给定初始值
xinitial = 1;
x = zeros(nely,nelx);
switch xinitial
    case 1 %中间一条1/5宽度的窄带
        len = floor(nely*volfrac); 
        lef = (nely-len)/2+1;
        rig = lef + len -1;
        fixeddofs = [lef : rig]';
        ind = repmat(fixeddofs,1,nelx)+repmat((0:nelx-1)*nely, len, 1);
        Ind = reshape(ind, len*nelx,1);
    case 2
        len = floor(sqrt(nely*nelx*volfrac));
        lef = (nely-len)/2+1;
        rig = lef + len -1;
        fixeddofs = [lef : rig]';
        ind = repmat(fixeddofs,1,len)+repmat(((lef-1):(rig-1))*nely, len, 1);
        Ind = reshape(ind, len^2,1);
end   
x(Ind) = 1;
if continuation == 1
    figure; imshow(1-x);
end
g = 0.7;
sd = 1;
[y,loop,c,x,energies]=topthr_penlity(nelx,nely,volfrac,lambda,g,sd,bc,continuation,x,fileID);
plot(1:loop, energies,'-')
xlabel('迭代次数');
ylabel('能量耗散');
hold on
figure; imshow(1-x); 