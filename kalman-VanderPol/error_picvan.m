clc; close all; clear all
t=0;
Tend=4999;
T=0.01;
for i=1:Tend   
    t(i+1)=(t(i)+T);
end
eva_error = dlmread('���Ƶ���error.txt');
act_error = dlmread('��������error.txt');
eva = dlmread('���Ƶ���.txt');
act = dlmread('��������.txt');
cegu=dlmread('��������error.txt');
% length_2= size(eva_error,1);
% plot(1:length_2,cegu(:,1),'r--')
% xlabel ('ts');ylabel ('\Deltax1/m');
% figure
% plot(1:length_2,cegu(:,2),'r--')
% xlabel ('ts');ylabel ('\Deltax2/m');
%  
    evax_error = eva_error(:,1); 
    evay_error = eva_error(:,2); 


    x_error = act_error(:,1);
    y_error = act_error(:,2);
    
    evax = eva(:,1);
    evay = eva(:,2);
    
    actx = act(:,1);
    acty = act(:,2);
%%%%%%%%%%%%%%%%%%%%%%%����������ʵ�ʵ����Ա���άͼ    
    plot(t,evax,'r-.');figure
    plot(t,actx,'b-.')
%     figure
%     plot(t,evay,'r-');
%     figure
%     plot(t,acty,'b-')
    %xlabel ('x');ylabel ('z');zlabel ('y');
%     grid on
%     figure(gcf + 1)
%     hold on
%     grid on
%     figure
%     plot(evax,'r.')
%     plot(actx,'r.')
   figure
    hold on

    plot(t,evay_error,'b--')
    figure
     plot(t,y_error+3,'b--')
   % xlabel ('x');ylabel ('y')
%%%%%%%%%%%%%%%%%%%%%%%�����귽���ϵ����Ա�ͼ   
