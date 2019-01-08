close all;clc;clear all;
global d
d = 98e-3;
global i43
i43 =4;
global mass
mass = 4.3;
global g
g = 9.8066;
global rho
rho = 1.2250;
global Sm
Sm = pi*d*d/4;
global Cx0
Cx0=0.157;
global T 
T=0.01;
global Cx
Cx = -rho*Cx0*i43/(2*mass)*Sm;
x(1)=0;y(1)=0;v(1)=110;theta(1)=25/180*pi;
t=0;
Tend=8.5/T;

f1 =  Cx*v*v-g*sin(theta);
f2 = - g*cos(theta)/v;
f3 = v*cos(theta);
f4 = v*sin(theta);

for i=1:Tend   
    v(i+1)= v(i)+T*(Cx*v(i)*v(i)-g*sin(theta(i)));
    theta(i+1)= theta(i)-T*g*cos(theta(i))/v(i);
    x(i+1)= x(i)+T*v(i)*cos(theta(i));
    y(i+1)= y(i)+T*v(i)*sin(theta(i));
    t(i+1)=(t(i)+T);
end

% F = [1 - Cx*Sm*T*mass*rho*v, -T*g*cos(theta), 0, 0;
%      (T*g*cos(theta))/v^2, (T*g*sin(theta))/v + 1, 0, 0;
%      T*cos(theta), -T*v*sin(theta), 1, 0;
%      T*sin(theta), T*v*cos(theta), 0, 1];
global H;
global G;
global R
H=eye(4);G=eye(4);
X=[v;theta;x;y];%每列表示一个状态向量数据组
R=[0.05,0,0,0;0,2/180*pi,0,0;0,0,5,0;0,0,0,5];
V1 = normrnd(0,sqrt(R(1,1)),1,Tend+1);
z1 = v+V1;
V2 = normrnd(0,sqrt(R(2,2)),1,Tend+1);
z2 = theta+V2;
V3 = normrnd(0,sqrt(R(3,3)),1,Tend+1);
z3 = x+V3;
V4 = normrnd(0,sqrt(R(4,4)),1,Tend+1);
z4 = y+V4;
z=[z1;z2;z3;z4];%每列表示一个测量数据组
% plot(z3,z4)
% figure
% plot(z3,z1);
% figure
% plot(z3,z2)
%p=n*n
% W=[0.0075,0,0,0;
%     0,0.0025,0,0;
%     0,0,0.0025,0
%     0,0,0,0.0045];
W=[0.0005,0,0,0;
    0,0.00005,0,0;
    0,0,0.5,0
    0,0,0,10];

Xest(:,1)=[110;20/180*pi;0;0];
p{1}=zeros(4);%4*4
dJ=0;
% X=[v;theta;x;y];%每列表示一个状态向量数据组
for i=1:Tend
    Xfest(:,i+1)=[Xest(1,i)+T*(Cx*Xest(1,i)*Xest(1,i)-g*sin(Xest(2,i)));
                  Xest(2,i)-T*g*cos(Xest(2,i))/Xest(1,i);
                  Xest(3,i)+T*Xest(1,i)*cos(Xest(2,i));
                  Xest(4,i)+T*Xest(1,i)*sin(Xest(2,i))];

    Mk{i+1}=(inv(R))*(z(:,i+1)-Xfest(:,i+1));
    Fy=[1 - Cx*Sm*T*mass*rho*Xfest(1,i+1), -T*g*cos(Xfest(2,i+1)), 0, 0;
       (T*g*cos(Xfest(2,i+1)))/Xfest(1,i+1)^2, (T*g*sin(Xfest(2,i+1)))/Xfest(1,i+1) + 1, 0, 0;
        T*cos(Xfest(2,i+1)), -T*Xfest(1,i+1)*sin(Xfest(2,i+1)), 1, 0;
        T*sin(Xfest(2,i+1)), T*Xfest(1,i+1)*cos(Xfest(2,i+1)), 0, 1];
    pk{i+1}=W+Fy*p{i}*Fy';
    p{i+1}=pk{i+1}*inv(eye(4)+Fy*pk{i+1});
    Xest(:,i+1)=Xest(:,i)+p{i+1}*Mk{i+1};
%     dJ=dJ+inv(eye(4)+Fy*pk{i+1})*Mk{i+1};
%     W=W+diag(0.05*dJ,0);
    STDVest(i+1)=std(z(1,1:i)-Xest(1,1:i));
    STDthest(i+1)=std(z(2,1:i)-Xest(2,1:i));
    STDXest(i+1)=std(z(3,1:i)-Xest(3,1:i));
    STDYest(i+1)=std(z(4,1:i)-Xest(4,1:i));
end
Xest;

% plot(x(1:820),y(1:820),'k')%true
% hold on
% plot(Xest(3,:),Xest(4,:),'r')%estimate
% hold on
% plot(z3(1:820),z4(1:820),'b')%measure
% legend('true','estimate','measure')
%  xlabel('射程/m');ylabel('射高/m')
% grid on
% plot(x,y,'k')%true
% hold on
% plot(Xest(3,:),Xest(4,:),'r')%estimate
% hold on
% plot(z3,z4,'b')%measure
% legend('true','estimate','measure')
% grid on
% 
% figure
% plot(x(1:820),v(1:820),'k',Xest(3,:),Xest(1,:),'r',z3(1:820),z1(1:820),'b')
% xlabel('射程/m');ylabel('速度/m/s')
% legend('true','estimate','measure')
% grid on
% figure

rV=R(1,1)*ones(1,Tend+1);
rth=R(2,2)*ones(1,Tend+1);
rX=R(3,3)*ones(1,Tend+1);
rY=R(4,4)*ones(1,Tend+1);
plot(t,STDVest.^2,t,rV)
xlabel('t');ylabel('速度方差')
axis([0 9 0 0.2])
figure
plot(t,STDthest.^2,t,rth)
xlabel('t');ylabel('角度方差')
axis([0 9 0 0.07])
figure
plot(t,STDXest.^2,t,rX)
xlabel('t');ylabel('x方差')
axis([0 9 0 70])
figure
plot(t,STDVest.^2,t,rY)
xlabel('t');ylabel('y方差')
% eV=Xest(1,:)-v;
% eth=Xest(2,:)-theta;
% ex=Xest(3,:)-x;
% ey=Xest(4,:)-y;
% mV=z(1,:)-v;
% mth=z(2,:)-theta;
% mx=z(3,:)-x;
% my=z(4,:)-y;
% plot(t,eV,'r',t,mV)
% xlabel('t');ylabel('速度偏差')
% figure
% plot(t,eth,'r',t,mth)
% xlabel('t');ylabel('角度偏差')
% figure
% plot(t,ex,'r',t,mx)
% xlabel('t');ylabel('x偏差')
% figure
% plot(t,ey,'r',t,my)
% xlabel('t');ylabel('y偏差')
% 
% % figure
% % plot(x(1:820),theta(1:820),'k',Xest(3,:),Xest(2,:),'r',z3(1:820),z2(1:820),'b')
% % xlabel('射程/m');ylabel('弹道倾角/rad')
% % legend('true','estimate','measure')
% % grid on
% % 
% % figure
% % plot(1:Tend+1,x,'k',1:Tend+1,Xest(3,:),'r',1:Tend+1,z3,'b')
% % legend('true','estimate','measure')
% % grid on
% % figure
% % plot(1:Tend+1,y,'k',1:Tend+1,Xest(4,:),'r',1:Tend+1,z4,'b')
% % legend('true','estimate','measure')
% % grid on
% % figure
% % plot(1:Tend+1,theta,'k',1:Tend+1,Xest(2,:),'r',1:Tend+1,z2,'b')
% % legend('true','estimate','measure')
% % grid on
% % figure
% % plot(1:Tend+1,v,'k',1:Tend+1,Xest(1,:),'r',1:Tend+1,z1,'b')
% % legend('true','estimate','measure')
% % grid on
% 
% 
% 
% 
% 
% 
% 
