% save the following output to this folder
save_to_folder = false;
% fileID for log file
fileID = 'penlity_general_method2.log';
logtype = 'a';
% whether calculate compliance with respect to the same Emin
sameEmin = 5e-5;
same = false;
% parameters
volfrac = 0.2;          %�����
nelx = 200;
nely = 200; 
Emin = [0.05,0.01,0.005,0.002];            %������ϵ��������ģ��
%bc = 'left_bdc_right_up_qin';
bc = 'left_down_bdc_right_down_qin';%�߽�ѡ
%bc = 'left_up_down_bdc_left_down_qin';
%objectfunc = 'down_central';
objectfunc = 'left_up';
%objectfunc =  'right_down';
w = 5; %�����Ȩ��
g = 0.0000;
sd = 1;           %sd/nely Ϊ�������tau
continuation = 0; %�Ƿ�ʹ��Ԥ�����״x
x = 1;            %Ԥ����״x
Vforce = 0.00;%�������Դ�С
lambda = 1;
r      = 10000; %�ڽ�����
%[y, loop, c,  x, energies, energies_k] = topthr_penlity(nelx, nely, lambda, volfrac, Emin(4), g(1), sd, objectfunc,bc, w,continuation, x, fileID,logtype,Vforce);
[y, loop, c,  x, energies, energies_k] = topthr_penlity_general(nelx, nely, lambda, r, volfrac, Emin(3), g, sd, objectfunc,bc, w,continuation, x, fileID,logtype,Vforce);
plot(1:loop, energies,'-')
xlabel('��������');
ylabel('Ŀ�꺯��');
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
%     xlabel('Ŀ�꺯��');
%     ylabel('������');
%     hold on
%     figure; imshow(1-x); 
% end