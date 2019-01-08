clc; close all; clear all
t=0;
Tend=4999;
T=0.01;
for i=1:Tend   
    t(i+1)=(t(i)+T);
end
eva_error = dlmread('估计弹道error.txt');
act_error = dlmread('测量弹道error.txt');
eva = dlmread('估计弹道.txt');
act = dlmread('测量弹道.txt');
cegu=dlmread('测量估计error.txt');
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
%%%%%%%%%%%%%%%%%%%%%%%测量弹道与实际弹道对比三维图    
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
%%%%%%%%%%%%%%%%%%%%%%%各坐标方向上的误差对比图   
