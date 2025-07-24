% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID = 'penlity_general_method2.log';
logtype = 'a';
% whether calculate compliance with respect to the same Emin
sameEmin = 5e-5;
same = false;
% parameters
volfrac = 0.2;          %体积比
nelx = 200;
nely = 200; 
Emin = [0.05,0.01,0.005,0.002];            %人造材料的相对杨氏模量
%bc = 'left_bdc_right_up_qin';
bc = 'left_down_bdc_right_down_qin';%边界选
%bc = 'left_up_down_bdc_left_down_qin';
%objectfunc = 'down_central';
objectfunc = 'left_up';
%objectfunc =  'right_down';
w = 5; %输出功权重
g = 0.0000;
sd = 1;           %sd/nely 为卷积参数tau
continuation = 0; %是否使用预设的形状x
x = 1;            %预设形状x
Vforce = 0.00;%体积力相对大小
lambda = 1;
r      = 10000; %邻近因子
%[y, loop, c,  x, energies, energies_k] = topthr_penlity(nelx, nely, lambda, volfrac, Emin(4), g(1), sd, objectfunc,bc, w,continuation, x, fileID,logtype,Vforce);
[y, loop, c,  x, energies, energies_k] = topthr_penlity_general(nelx, nely, lambda, r, volfrac, Emin(3), g, sd, objectfunc,bc, w,continuation, x, fileID,logtype,Vforce);
plot(1:loop, energies,'-')
xlabel('迭代次数');
ylabel('目标函数');
hold on
figure; imshow(1-x); 
% for i = 1:length(Emin)
%     if i > 1
%         continuation = 1;
%         logtype = 'a';
%     end
%     %[y, loop, c,  x, energies, energies_k] = topthr_penlity(nelx, nely, lambda, volfrac, Emin(i), g(1), sd, objectfunc,bc, w,continuation, x, fileID,logtype,Vforce);
%     [y, loop, c,  x, energies, energies_k] = topthr_penlity_general(nelx, nely, lambda, volfrac, Emin(4), g(1), sd, objectfunc,bc, w,continuation, x, fileID,logtype,Vforce);
%     plot(1:loop, energies,'-')
%     xlabel('目标函数');
%     ylabel('柔性能');
%     hold on
%     figure; imshow(1-x); 
% end