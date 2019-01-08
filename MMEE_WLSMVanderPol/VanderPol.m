 close all;
clear all;
clc;
x1(1)=-1;
x2(1)=-1;
T=0.01;
t=0;
Tend=4999;
for i=1:Tend   
    x1(i+1)=(1+T)*x1(i)-T*(x2(i)+x1(i)*x2(i)*x2(i));
    x2(i+1)=x2(i)+T*x1(i);
    t(i+1)=(t(i)+T);
end

b=0;
r=1;
v = normrnd(0,sqrt(r),1,Tend+1);
z=x1+b+v;
 
w=[0.005 0;0 0.05];
Xest(:,1)=[0.35;0.35];
% Xest(:,1)=[0.25;0.25];
P{1}=[0,0;0,0];
dJ=0;
% Xest(:,1)=[X(1,1)-T*(X(2,1)+X(1,1)*X(2,1)*X(2,1));
%        X(2,1)+T*X(1,1)];
% M(1)=1/r*(z(1)-Xest(1,1));
% F=[1-T*Xest(2,1) -T-2*T*Xest(1,1)*Xest(2,1);T 1];
% Pk{1}=w+F*P{1}*F';
% P{1}=Pk{1}*inv(eye(2)+F*Pk{1})
% Xest(:,2)=Xest(:,1)+P{1}*[M(1);0];
for i=1:Tend
    Xfest(:,i+1)=[Xest(1,i)-T*(Xest(2,i)+Xest(1,i)*Xest(2,i)*Xest(2,i));
               Xest(2,i)+T*Xest(1,i)];
    M(i+1)=1/r*(z(i+1)-Xfest(1,i+1));
    F=[1-T*Xfest(2,i+1) -T-2*T*Xfest(1,i+1)*Xfest(2,i+1);
        T 1];
    Pk{i+1}=w+F*P{i}*F';
    P{i+1}=Pk{i+1}*inv(eye(2)+F*Pk{i+1});
    Xest(:,i+1)=Xest(:,i)+P{i+1}*[M(i+1);0];
    STDest(i+1)=std(z(1:i)-Xest(1,1:i));
%     dJ= inv(eye(2)+F*Pk{i+1})*[M(i+1)];
%     if (r-STDest(i+1)^2)>0.1
%     w=w-0.01*dJ;
%     end
end
 
plot(t,z,'g')%%测量值 绿色
% hold on
figure
plot(t,x1)%%真实值 蓝色
hold on
% figure
plot(t,Xest(1,:),'r')%%估计值，红色
figure
R=r*ones(1,Tend+1);%%测量误差方差
plot(t,STDest,t,R)
figure
plot(t,Xest(1,:)-x1,'r')%%估计值，红色
figure
plot(t,x2)
hold on
plot(t,Xest(2,:),'r')








