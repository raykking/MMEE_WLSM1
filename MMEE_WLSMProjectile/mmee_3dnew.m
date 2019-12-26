function mmee_3dnew
clear;close all;clc;
digits(5) 
true_data = dlmread('true_values.txt');
length_1 = size(true_data,1);
% for j = 1:length_1
%     data(j,5) = data(j,5)/180*pi;
%     data(j,6) = data(j,6)/180*pi;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms x y z v theta phi nx ny nz T;
% % nx ny nz are linear error of x,y,z that we want to estimate; a trial
% Wx = 0;%%%%%%%%%%%%%%%%%%%%%%%%%
% Wz = 0;%%%%%%%%%%%%%%%%%%%%%%%
% g = 9.8066;
% Cx0 = 0.157;rho = 1.2250;i43 =4;d = 98e-3;mass=4.3;
% Sm = pi*d*d/4;
% Cx = -rho*Cx0*i43/(2*mass)*Sm;
% wx = Wx*cos(phi)*cos(theta)+Wz*sin(phi);
% wy = -Wx*sin(theta);
% wz = -Wx*sin(phi)*cos(theta)+Wz*cos(phi);
% vrr = (v-wx)*(v-wx)+wy*wy+wz*wz;
% vr = vrr^(1/2); 

%%%%Yakebi_ matrix
% f(x,y,z,v,theta,phi,nx,ny,nz)=[v*cos(theta)*cos(phi);
%                                 v*sin(theta)*cos(phi);
%                                 v*sin(phi);
%                                 Cx*vr*(v-wx)-g*sin(theta)*cos(phi);
%                                 (-Cx*vr*wy-g*cos(theta))/(v*cos(phi));
%                                 (-Cx*vr*wz+g*sin(theta)*sin(phi))/v;
%                                 0;
%                                 0;
%                                 0;];
% J(x,y,z,v,theta,phi,nx,ny,nz) = jacobian(...
%                                 [v*cos(theta)*cos(phi);
%                                  v*sin(theta)*cos(phi);
%                                  v*sin(phi);
%                                  Cx*vr*(v-wx)-g*sin(theta)*cos(phi);
%                                  (-Cx*vr*wy-g*cos(theta))/(v*cos(phi));
%                                  (-Cx*vr*wz+g*sin(theta)*sin(phi))/v;
%                                  0;
%                                  0;
%                                  0],[x y z v theta phi nx ny nz]);
%%%% J = df/dX
% J(x,y,z,v,theta,phi,nx,ny,nz) = jacobian(...
%                                 [v*cos(theta)*cos(phi);
%                                  v*sin(theta)*cos(phi);
%                                  v*sin(phi);
%                                  Cx*vr*(v-(Wx*cos(phi)*cos(theta)+Wz*sin(phi)))-g*sin(theta)*cos(phi);
%                                  (-Cx*vr*(-Wx*sin(theta))-g*cos(theta))/(v*cos(phi));
%                                  (-Cx*vr*(-Wx*sin(phi)*cos(theta)+Wz*cos(phi))+g*sin(theta)*sin(phi))/v;
%                                  0;
%                                  0;
%                                  0],[x y z v theta phi nx ny nz]);
% f(x,y,z,v,theta,phi,nx,ny,nz) =[v*cos(theta)*cos(phi);
%     v*sin(theta)*cos(phi);
%     v*sin(phi);
%     Cx*v*(v-0)-g*sin(theta)*cos(phi);
%     (-Cx*v*0-g*cos(theta))/(v*cos(phi));
%     (-Cx*v*0+g*sin(theta)*sin(phi))/v;
%     0;
%     0;
%     0];

% J(x,y,z,v,theta,phi,nx,ny,nz) = vpa(simplify(jacobian(...
%                                 [v*cos(theta)*cos(phi);
%                                  v*sin(theta)*cos(phi);
%                                  v*sin(phi);
%                                  Cx*v*(v-0)-g*sin(theta)*cos(phi);
%                                  (-Cx*v*0-g*cos(theta))/(v*cos(phi));
%                                  (-Cx*v*0+g*sin(theta)*sin(phi))/v;
%                                  0;
%                                  0;
%                                  0],[x y z v theta phi nx ny nz])),3);
%                              
%                              
                    
                             
% J(x,y,z,v,theta,phi,nx,ny,nz) = [ 0, 0, 0,                cos(phi)*cos(theta),       -1.0*v*cos(phi)*sin(theta),                   -1.0*v*cos(theta)*sin(phi), 0, 0, 0;
%   0, 0, 0,                cos(phi)*sin(theta),            v*cos(phi)*cos(theta),                   -1.0*v*sin(phi)*sin(theta), 0, 0, 0;
%   0, 0, 0,                           sin(phi),                                0,                                   v*cos(phi), 0, 0, 0;
%   0, 0, 0,                       -0.0013495*v,      -9.8066*cos(phi)*cos(theta),                   9.8066*sin(phi)*sin(theta), 0, 0, 0;
%   0, 0, 0, (9.8066*cos(theta))/(v^2*cos(phi)), (9.8066*sin(theta))/(v*cos(phi)), -(9.8066*cos(theta)*sin(phi))/(v*cos(phi)^2), 0, 0, 0;
%   0, 0, 0,  -(9.8066*sin(phi)*sin(theta))/v^2,   (9.8066*cos(theta)*sin(phi))/v,               (9.8066*cos(phi)*sin(theta))/v, 0, 0, 0;
%   0, 0, 0,                                  0,                                0,                                            0, 0, 0, 0;
%   0, 0, 0,                                  0,                                0,                                            0, 0, 0, 0;
%   0, 0, 0,                                  0,                                0,                                            0, 0, 0, 0];
% 
% f(x,y,z,v,theta,phi,nx,ny,nz) = [x;y;z;v;theta;phi;nx;ny;nz] + T * f(x,y,z,v,theta,phi,nx,ny,nz) + 0.5 * T^2 * J(x,y,z,v,theta,phi,nx,ny,nz) * f(x,y,z,v,theta,phi,nx,ny,nz);
% surr_f(x,y,z,v,theta,phi,nx,ny,nz) = vpa(simplify(f(x,y,z,v,theta,phi,nx,ny,nz)),3)
% dsurr_f = jacobian(surr_f,[x y z v theta phi nx ny nz])

% Jf_f = vpa(simplify(Jf*f),3); %9*1,df/dX * f
% G(k)
% G = [zeros(6,6);ones(3,6)];
G = eye(6,6);
% assume d = sqrt(D)*randn(3,1) d is what we want to estimate
 
% linear system model
% ff(xk) = xk + T*f(xk)+0.5*T^2*J*f(xk)

% X(k+1) = ff(x(k))+G(k)*d(k)
% Measurement noise=Rc,sigma^2 of x,y,z,v,theta, phi
% Z(k) = h(X(k))+e, the variance of e is R
% we set h(x) is linear, equals to I9*9
R = [3 0 0 0 0 0;
    0 6 0 0 0 0;
    0 0 3 0 0 0;
    0 0 0 5 0 0;
    0 0 0 0 0.01 0;
    0 0 0 0 0 0.01];
% h(x,y,z,v,theta,phi,nx,ny,nz) = [x;y;z;v;theta;phi];
% h 6*1
% H(x,y,z,v,theta,phi,nx,ny,nz) = jacobian([x;y;z;v;theta;phi],[x,y,z,v,theta,phi,nx,ny,nz])
H = eye(6);
%H 6*6
x = true_data(:,1);
y = true_data(:,2);
z = true_data(:,3);
v = true_data(:,4);
theta = true_data(:,5);
phi = true_data(:,6);
 
% H=eye(6);G=eye(6);
X=[x,y,z,v,theta,phi];% true value
V1 = normrnd(0,sqrt(R(1,1)),length_1,1);
z1 = x+V1;

V2 = normrnd(0,sqrt(R(2,2)),length_1,1);
z2 = y+V2;

V3 = normrnd(0,sqrt(R(3,3)),length_1,1);
z3 = z+V3;

V4 = normrnd(0,sqrt(R(4,4)),length_1,1);
z4 = v+V4;

V5 = normrnd(0,sqrt(R(5,5)),length_1,1);
z5 = theta+V5;

V6 = normrnd(0,sqrt(R(6,6)),length_1,1);
z6 = phi+V6;
Z=[z1,z2,z3,z4,z5,z6];% N*6

% uodate:
% F = jacobian(ff(X))
% F(x,y,z,v,theta,phi,nx,ny,nz) = jacobian(,[(x,y,z,v,theta,phi,nx,ny,nz)])
% M(k+1)=(H(k+1)'/(R))*(Z(k+1)-h(Xest_k(k+1)))
% Xest_k(k+1) = ff(X_est(k))
% Pk(k+1)=G(k)*W*G(k)'+F(Xest_k(k+1))*P(k)*F(Xest_k(k+1))'
% Mx = jacobian(M(X))
% P(k+1)=Pk(k+1)/((I-Mx(Xest_k(k+1))*Pk(k+1)))
% Xest(k+1) = ff(X_est(k))+P(k+1)*M(k+1);
% Weight matrix: W
W = [0.00008,0,0,0,0,0;
     0,0.0000008,0,0,0,0;
     0,0, 0.0000008,0,0,0;
     0,0,0, 0.0000008,0,0;
     0,0,0,0, 0.0000008,0;
     0,0,0,0,0,0.0000008];
X_est(:,1)=[x(1),y(1),z(1)-1,v(1),theta(1),phi(1)]';
% X_est(:,1) = [0,0,0,120,0.7,0];
T = 0.01;
I9 = eye(6);
P{1} = zeros(6,6);
f_X_est(:,1) = surr_f(X_est(:,1), T);
X_est_k = f_X_est(:,1); 
M{1} = H' / R * (Z(1,:)' - h(X_est_k));
dM = - surr_F(X_est_k,T);
Pk{1} = G * W * G' + surr_F(X_est_k,T) * P{1} * surr_F(X_est_k,T)';
P{1} = Pk{1} / (I9 - dM * Pk{1});


for k = 1:length_1-1
    f_X_est(:,k+1) = surr_f(X_est(:,k), T);
    
    % X_est(k+1|k)
    X_est_k = f_X_est(:,k+1); 
    
    % M(k+1)
    M{k+1} = H' / R * (Z(k+1,:)' - h(X_est_k));  % 6*1
    dM = - surr_F(X_est_k,T);
    % P(k+1|k)
    Pk{k+1} = G * W * G' + surr_F(X_est_k,T) * P{k} * surr_F(X_est_k,T)';
    
    % P(k+1)
    P{k+1} = Pk{k+1} / (I9 - dM * Pk{k+1});
  
    % X estimate next
    X_est(:,k+1) = f_X_est(:,k+1) + P{k+1} * M{k+1};
    Xest = X_est';
    STDVest(k+1) = std(Z(1:k,4) - Xest(1:k,4));
    STDthest(k+1) = std(Z(1:k,5) - Xest(1:k,5));
    STDphiest(k+1) = std(Z(1:k,6) - Xest(1:k,6));
    STDXest(k+1) = std(Z(1:k,1) - Xest(1:k,1));
    STDYest(k+1) = std(Z(1:k,2) - Xest(1:k,2));
    STDZest(k+1) = std(Z(1:k,3) - Xest(1:k,3));
end

fp = fopen('est.txt','w');
save('est.txt','Xest','-ascii','-append');
fp = fopen('celi.txt','w');
save('celi.txt','Z','-ascii','-append');
% fp = fopen('true-est-error.txt','w');
% save('true-est-error.txt','data_1','-ascii','-append');
% fp = fopen('true-celi-error.txt','w');
% save('true-celi-error.txt','delta_act1','-ascii','-append');
plot3(x,z,y,'k','LineWidth',1)
hold on 

plot3(Xest(:,1),Xest(:,3),Xest(:,2),'r--','LineWidth',1)
hold on

plot3(z1,z3,z2,'b','LineWidth',0.5)
grid on
xlabel('{\itx}/m');ylabel('{\itz}/m');zlabel('{\ity}/m')
legend('true','estimate','measurement')
figure
rX=R(1,1)*ones(1,length_1);
rY=R(2,2)*ones(1,length_1);
rZ=R(3,3)*ones(1,length_1);
rV=R(4,4)*ones(1,length_1);
rth=R(5,5)*ones(1,length_1);
rPhi=R(6,6)*ones(1,length_1);
t = 1:length_1;


subplot(3,2,1);
plot(t/100,STDXest.^2,'r--',t/100,rX,'b','LineWidth',1)
xlabel('{\itt}/s');ylabel('{\itrange variance}/m')
grid on
% figure
subplot(3,2,2);
plot(t/100,STDYest.^2,'r--',t/100,rY,'b','LineWidth',1)
xlabel('{\itt}/s');ylabel('{\itheight variance}/m')
grid on
% figure
subplot(3,2,3);
plot(t/100,STDZest.^2,'r--',t/100,rZ,'b','LineWidth',1)
xlabel('{\itt}/s');ylabel('{\itlateral range variance}/m')
grid on
% figure
subplot(3,2,4);
plot(t/100,STDVest.^2,'r--',t/100,rV,'b','LineWidth',1)
xlabel('{\itt}/s');ylabel('{\itvelocity variance}/m/s')
grid on 
% figure
subplot(3,2,5);
plot(t/100,STDthest.^2,'r--',t/100,rth,'b','LineWidth',1)
xlabel('{\itt}/s');ylabel('{\it\theta variance}/rad')
grid on 
% figure
subplot(3,2,6);
plot(t/100,STDphiest.^2,'r--',t/100,rPhi,'b','LineWidth',1)
xlabel('{\itt}/s');ylabel('{\it\phi variance}/rad')
grid on
figure
eV=Xest(:,4)-v;
eth=Xest(:,5)-theta;
ephi=Xest(:,6)-phi;
ex=Xest(:,1)-x;
ey=Xest(:,2)-y;
ez=Xest(:,3)-z;

mV=Z(:,4)-v;
mth=Z(:,5)-theta;
mphi=Z(:,6)-phi;
mx=Z(:,1)-x;
my=Z(:,2)-y;
mz=Z(:,3)-z;

% figure
% plot((1:850)/100,ex(1:850),'r--')
% xlabel('{\itt}/s');ylabel('x-error');grid on 

subplot(3,2,1);
plot(t/100,ex,'r--',t/100,mx,'b')
xlabel('{\itt}/s');ylabel('x-error');grid on 
subplot(3,2,2);
plot(t/100,ey,'r--',t/100,my,'b')
xlabel('{\itt}/s');ylabel('y-error');grid on 
subplot(3,2,3);
plot(t/100,ez,'r--',t/100,mz,'b')
xlabel('{\itt}/s');ylabel('z-error');grid on 
subplot(3,2,4);
plot(t/100,eV,'r--',t/100,mV,'b')
xlabel('{\itt}/s');ylabel('velocity-error');grid on 
subplot(3,2,5);
plot(t/100,eth,'r--',t/100,mth,'b')
xlabel('{\itt}/s');ylabel('theta-error');grid on 
subplot(3,2,6);
plot(t/100,ephi,'r--',t/100,mphi,'b')
xlabel('{\itt}/s');ylabel('phi-error');grid on 

legend('estimate-error','measurement-error')

 

function y = h(X)
    y = [X(1);X(2);X(3);X(4);X(5);X(6)];  
% f
function y = surr_f(X,T)
    g = 9.8066;
    Cx0 = 0.157;rho = 1.2250;i43 =4;d = 98e-3;mass=4.3;
    Sm = pi*d*d/4;
    Cx = -rho*Cx0*i43/(2*mass)*Sm;
 
    y = [                                                                                                                                                        - 0.00033737*cos(X(6))*cos(X(5))*T^2*X(4)^2 + cos(X(6))*cos(X(5))*T*X(4) + X(1);
                                                                                                                                                    - 0.00033737*cos(X(6))*sin(X(5))*T^2*X(4)^2 - 4.9033*T^2 + cos(X(6))*sin(X(5))*T*X(4) + X(2);
                                                                                                                                                                                       - 0.00033737*sin(X(6))*T^2*X(4)^2 + sin(X(6))*T*X(4) + X(3);
                                        X(4) - 1.0*T*(0.00067474*X(4)^2 + 9.8066*cos(X(6))*sin(X(5))) + (48.085*T^2*cos(X(5))^2)/X(4) + 0.00067474*T^2*X(4)*(0.00067474*X(4)^2 + 9.8066*cos(X(6))*sin(X(5))) + (48.085*T^2*sin(X(6))^2*sin(X(5))^2)/X(4);
 -(3.4694e-26*(2.7719e+27*T^2*cos(X(5))*sin(X(5)) - 2.8823e+25*X(5)*X(4)^2*cos(X(6))^2 + 5.1305e+19*T^2*cos(X(6))^2*cos(X(5))*sin(X(5)) + 9.536e+22*T^2*X(4)^2*cos(X(6))*cos(X(5)) + 2.8266e+26*T*X(4)*cos(X(6))*cos(X(5))))/(X(4)^2*cos(X(6))^2);
  (1.0e-8*(1.0e+8*X(6)*X(4)^2*cos(X(6)) - 4.8085e+9*T^2*cos(X(5))^2*sin(X(6)) + 9.6169e+9*T^2*cos(X(6))^2*sin(X(6))*sin(X(5))^2 + 330850.0*T^2*X(4)^2*cos(X(6))*sin(X(6))*sin(X(5)) + 9.8066e+8*T*X(4)*cos(X(6))*sin(X(6))*sin(X(5))))/(X(4)^2*cos(X(6)))];
 

% df/dx                                                                                                                                                                         
function y = surr_F(X,T)
    g = 9.8066;
    Cx0 = 0.157;rho = 1.2250;i43 =4;d = 98e-3;mass=4.3;
    Sm = pi*d*d/4;
    Cx = -rho*Cx0*i43/(2*mass)*Sm;
    y = [ 1.0,   0,   0,                                                                                                         -1.08e-19*T*cos(X(6))*cos(X(5))*(6.22e+15*T*X(4) - 9.22e+18),                                                                                                                                                                                5.42e-20*T*X(4)*cos(X(6))*sin(X(5))*(6.22e+15*T*X(4) - 1.84e+19),                                                                                                                                                                                                                            5.42e-20*T*X(4)*cos(X(5))*sin(X(6))*(6.22e+15*T*X(4) - 1.84e+19);
            0, 1.0,   0,                                                                                                         -1.08e-19*T*cos(X(6))*sin(X(5))*(6.22e+15*T*X(4) - 9.22e+18),                                                                                                                                                                               -5.42e-20*T*X(4)*cos(X(6))*cos(X(5))*(6.22e+15*T*X(4) - 1.84e+19),                                                                                                                                                                                                                            5.42e-20*T*X(4)*sin(X(6))*sin(X(5))*(6.22e+15*T*X(4) - 1.84e+19);
            0,   0, 1.0,                                                                                                                    -1.08e-19*T*sin(X(6))*(6.22e+15*T*X(4) - 9.22e+18),                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                      -5.42e-20*T*X(4)*cos(X(6))*(6.22e+15*T*X(4) - 1.84e+19);
            0,   0,   0, 6.75e-4*T^2*(6.75e-4*X(4)^2 + 9.81*cos(X(6))*sin(X(5))) - 0.00135*T*X(4) + 9.11e-7*T^2*X(4)^2 - (48.1*T^2*cos(X(5))^2)/X(4)^2 - (48.1*T^2*sin(X(6))^2*sin(X(5))^2)/X(4)^2 + 1.0,                                                                                                                                      -(2.17e-23*T*cos(X(6))*cos(X(5))*(- 3.05e+20*T*X(4)^2 + 4.52e+23*X(4) + 4.44e+24*T*cos(X(6))*sin(X(5))))/X(4),                                                                                                                                                                                   (2.17e-23*T*sin(X(6))*sin(X(5))*(- 3.05e+20*T*X(4)^2 + 4.52e+23*X(4) + 4.44e+24*T*cos(X(6))*sin(X(5))))/X(4);
            0,   0,   0,                                         (2.27e-21*T*cos(X(5))*(8.46e+22*T*sin(X(5)) + 4.31e+21*X(4)*cos(X(6)) + 1.57e+15*T*cos(X(6))^2*sin(X(5))))/(X(4)^3*cos(X(6))^2), (7.52e-37*(2.37e+30*T^2*cos(X(6))^2 - 2.56e+38*T^2*cos(X(5))^2 + 1.33e+36*X(4)^2*cos(X(6))^2 + 1.28e+38*T^2 - 4.73e+30*T^2*cos(X(6))^2*cos(X(5))^2 + 4.4e+33*T^2*X(4)^2*cos(X(6))*sin(X(5)) + 1.3e+37*T*X(4)*cos(X(6))*sin(X(5))))/(X(4)^2*cos(X(6))^2),                                                                                                                                                            -(1.16e-18*T*cos(X(5))*sin(X(6))*(2.84e+15*T*cos(X(6))*X(4)^2 + 8.42e+18*cos(X(6))*X(4) + 1.65e+20*T*sin(X(5))))/(X(4)^2*cos(X(6))^3);
            0,   0,   0,            -(2.0e-4*T*sin(X(6))*(9.62e+5*T*cos(X(6))^2 - 4.81e+5*T*cos(X(5))^2 - 9.62e+5*T*cos(X(6))^2*cos(X(5))^2 + 4.9e+4*X(4)*cos(X(6))*sin(X(5))))/(X(4)^3*cos(X(6))),                                                                                       (5.0e-7*T*cos(X(5))*sin(X(6))*(1.92e+8*T*sin(X(5)) + 1.96e+7*X(4)*cos(X(6)) + 6620.0*T*X(4)^2*cos(X(6)) + 3.85e+8*T*cos(X(6))^2*sin(X(5))))/(X(4)^2*cos(X(6))), (5.0e-7*(3.85e+8*T^2*cos(X(6))^4 - 1.92e+8*T^2*cos(X(6))^2 - 9.62e+7*T^2*cos(X(5))^2 + 2.0e+6*X(4)^2*cos(X(6))^2 + 1.92e+8*T^2*cos(X(6))^2*cos(X(5))^2 - 3.85e+8*T^2*cos(X(6))^4*cos(X(5))^2 + 6620.0*T^2*X(4)^2*cos(X(6))^3*sin(X(5)) + 1.96e+7*T*X(4)*cos(X(6))^3*sin(X(5))))/(X(4)^2*cos(X(6))^2)];
    
 






