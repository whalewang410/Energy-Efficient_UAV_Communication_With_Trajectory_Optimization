clear all;
close all;
clc;
%Main_Constrained_trajectory
%1.UAV 参数
H = 100; %m
c1 = 9.26*1e-4;
c2 = 2250;
q0 = [0,1000].';%初始位置
qF = [1000,0].';%结束位置
v_0F = (qF-q0)/norm(qF-q0,2); %%初始和结束的单位速度方向
v0 = 30*v_0F;   %初始速度 m/s
vF = v0;        %结束速度
Vmax = 100;     %最大速度
Vmin = 3;       %最小速度，保证在空中飞行
amax = 5;      %最大加速度 m/s^2
T = 400;        %观察周期
deltat = 0.2;   %离散采样间隔
g = 9.8;        %重力加速度
N = T/deltat;
%2.Communication 参数
B = 1e6;        %带宽 MHz
N0dBm = -170;    %噪声功率谱dBm/Hz
N0 = 10^(N0dBm/10);
sigma2 = N0*B;  %噪声功率
PdBm = 10;      %无人机发射功率（固定）
P = 10^(PdBm/10);
beta0dB = -50;  %1m参考距离对应路损
beta0 = 10^(beta0dB/10);
gama0 = P*beta0/sigma2;

%3.Trajectory optimization
[q_EE,v_EE,a_EE,R_aveEE,P_aveEE,EE] = EE_A1(B,gama0,H,g,c1,c2,v0,vF,q0,qF,Vmax,Vmin,amax,deltat,T);
[q_RM,v_RM,a_RM,R_aveRM,P_aveRM,EE_RM] = RM_A1(B,gama0,H,g,c1,c2,q_EE,v_EE,v0,vF,q0,qF,Vmax,Vmin,amax,deltat,T);
[q_EM,v_EM,a_EM,R_aveEM,P_aveEM,EE_EM] = ...
    EM_A1(B,gama0,H,g,c1,c2,q_EE,v_EE,v0,vF,q0,qF,Vmax,Vmin,amax,deltat,T);

%Trajectory
figure(2)
plot(q_EE(1,:),q_EE(2,:),'k-','LineWidth',1.5)
hold on,grid on;
plot(q_RM(1,:),q_RM(2,:),'b-.','LineWidth',1.5)
plot(q_EM(1,:),q_EM(2,:),'g--','LineWidth',1.5)
title('Trajectory')
legend('EE-max','Rate-max','Energy-min')
xlabel('x(m)')
ylabel('y(m)')

%speed v.s. t
figure(3)
plot((1:N+1)*deltat,norms(v_EE,2),'k-','LineWidth',1.5);
hold on,grid on;
plot((1:N+1)*deltat,norms(v_RM,2),'b-.','LineWidth',1.5);
plot((1:N+1)*deltat,norms(v_EM,2),'g--','LineWidth',1.5);
title('Speed')
legend('EE-max','Rate-max','Energy-min')
xlabel('time(s)')
ylabel('Speed(m/s)')
xlim([0,200])
