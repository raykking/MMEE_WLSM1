VanderPol
%vanderPol
%x1=x';
%x2=x;
%%%x1'=(1-x2^2)*x1-x2;
%%%x2'=x1
clear all;close all;clc;
x1(1)=0.35;
x2(1)=0.35;
T=0.01;
t=0;
Tend=4999;
for i=1:Tend   
    x1(i+1)=(1+T)*x1(i)-T*(x2(i)+x1(i)*x2(i)*x2(i));
    x2(i+1)=x2(i)+T*x1(i);
    t(i+1)=(t(i)+T);
end
data=[x1' x2'];
length_1 = max(size(x1));
ts = 0.01;
xx1 = [];
xx2 = [];

delta_act = [];
X_act=[];
%%%%Process noise
Q=[1 0;0 1];
% W = sqrt(Q)*randn(6,1);
%%%%Measurement noise
R = [1 0;0 1];
% d = sqrt(R)*randn(6,1);
P = [100 0;0 100];
    
I1 = eye(2);
I2 = eye(2,2);
syms x1 x2
% FF(x1,x2)=[(1-x2^2)*x1-x2;x1];%û��ϵͳ���
FF(x1,x2)=[(-x2^2)*x1-x2;x1];%����ϵͳ���
% J(x1,x2) = jacobian([(1-x2^2)*x1-x2;x1],[x1 x2]);
J(x1,x2) = jacobian([(-x2^2)*x1-x2;x1],[x1 x2]);
 


for i = 1:length_1
    W = sqrt(Q)*randn(2,1);%ģ�����
    d = sqrt(R)*randn(2,1);%�۲����
%     ST=data(k,:)';%X(k-1)ʱ������ֵ
    
    x1 = data(i,1);
    x2 = data(i,2);
 
    f = J(x1,x2);
    ff = FF(x1,x2)+W;
    F = vpa(f,3);
    X = data(i,:)';%%%%%%��ʵ����
    %     X = a;%%%%%%��ʵ����
    h = data(i,:)';
    Z = h + d+3;%%%%%%%%����ֵ
    delta_act(:,end+1)=Z;
    %%%%
    Fai = I1 + F * ts;
    %%%%%%%%%״̬Ԥ��
    %     X = Fai*X;
    X = X + ff*ts;
    h = X;
    %%%%%%%%%%%%%Э����Ԥ��
    P = Fai*P*Fai'+ Q*ts;
    H = eye(2,2);
    %%%%%%%%%%RΪR(k+1),KΪ����
    Inverse_1 = H*P*H'+ R;
    Inverse_ = inv(Inverse_1);
    K = P*H'*Inverse_;
    %%%%%״̬����
    X = X+K*(Z-h);
    % %%%%%%%Э�������
    %    XX(end+1,1:5) = X';
    P = (I2 - K*H)*P;
    
    xx1(end+1,:) =X(1);
    xx2(end+1,:) =X(2);

 end
plot(t,xx2)
%
    data_ = [xx1,xx2];%%%%%%%%%���Ƶ�������data_ 
 
   delta_act = delta_act';%%%%%%%%������������delta_act
        for k = 1:2
            delta_act1(:,k) = data(:,k)- delta_act(:,k);%%%%%%%%%%%%%%%�����������delta_act1
            data_1(:,k) = data(:,k)-data_(:,k);%%%%%%%%%���Ƶ������data_1
            delta_cegu(:,k)=delta_act(:,k)-data_(:,k);
        end
        
 
   
       
fp = fopen('���Ƶ���.txt','w');
save('���Ƶ���.txt','data_','-ascii','-append');
fp = fopen('��������.txt','w');
save('��������.txt','delta_act','-ascii','-append');
fp = fopen('���Ƶ���error.txt','w');
save('���Ƶ���error.txt','data_1','-ascii','-append');
fp = fopen('��������error.txt','w');
save('��������error.txt','delta_act1','-ascii','-append');
fp = fopen('��������error.txt','w');
save('��������error.txt','delta_cegu','-ascii','-append');
%  error_picvan  