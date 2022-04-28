function [q,v,a,Raver,Paver,EE] = EE_A1(B,gama0,H,g,c1,c2,v0,vF,q0,qF,Vmax,Vmin,amax,deltat,T)

%Trajectory optimization to achieve energy-efficiency maximization design

%Input
%   B:         Bandwidth
%   gama0:     Reference received signal-to-noise ratio (SNR) at d0 = 1m.
%   H:         Flight Height
%   g:         Gravitational acceleration
%   c1:        UAV parameter 
%   c2:        
%   v0:        Initial velocity
%   vF:        Final velocity
%   q0:        Initial location
%   qF:        Final location
%   Vmax:
%   Vmin:
%   amax:      Max acceleration
%   deltat:    Time discretization interval
%   T:         Operation time
%Output:
%   q:          Location for different slot
%   v:          Velocity for different slot
%   a:          Acceleration for dofferent slot
%   Raver:      Average rate
%   Paver:      Average power
%   EE:         Energy efficiency

%1.Initialize 
%Straight flight with const speed as initial trajectory
N = T/deltat;   %Total slots
d = qF-q0;
v_0 = zeros(2,N+1); %N slots,N+1 speeds,x and y direction
q_0 = zeros(2,N+1);
a_0 = zeros(2,N);
v_aver = d/T;
v_0 = v_0+v_aver;
q_0(:,1) = q0;
for n = 2:N+1
    q_0(:,n) = q_0(:,n-1)+v_0(:,n-1)*deltat;
end
%——————————————————————————————————
%uniformly accelerated and retarded rectilinear motion 不能作为初始化点
% N = T/deltat;   %Total slots
% d = qF-q0;
% Tmid = T/2;
% vconst = (d-v0*deltat)./(T-3*deltat);
% a0 = (vconst-v0)/deltat;
% dmid = d/2;
% vmid = dmid/Tmid*2-v0;
% amid = (vmid-v0)/Tmid;
% a_0(:,1:N/2) = ones(2,N/2).*amid;    %Acceleration of first half path
% a_0(:,N/2+1:N) = -ones(2,N/2).*amid;
% 
% for n = 1:N+1
%     if n == 1
%         v_0(:,n) = v0;
%         q_0(:,n) = q0;
%         continue;
%     end
%     v_0(:,n) = v_0(:,n-1)+a_0(:,n-1)*deltat;
%     q_0(:,n) = q_0(:,n-1)+v_0(:,n-1)*deltat+1/2*a_0(:,n-1)*deltat^2;
% end
%——————————————————————————————————————
obj_SCA = zeros(1000,1);      %经过SCA后的目标函数值随迭代变化
obj_eUB = zeros(1000,1);      %能耗上界的目标函数值随迭代变化
obj_ee = zeros(1000,1);       %精确目标函数值随迭代变化
epsilonEE = 1;
th = 1e-3;
%2.SCA and Quadratic Transform
q_j = q_0;
v_j = v_0;
a_j = a_0;
t = 0;
obj1 = 0;
while epsilonEE >= th
    %SCA
    t = t+1;
    alpha_j = log2(1+gama0./(H^2+vecnorm(q_j(:,1:N),2).^2)).';
    beta_j = (gama0./(log(2)*(H^2+gama0+vecnorm(q_j(:,1:N),2).^2).*(H^2+vecnorm(q_j(:,1:N),2).^2))).';
    %Quadratic transform
    [q,v,a,obj] = QT(g,B,alpha_j,beta_j,c1,c2,v0,vF,q0,qF,v_j,q_j,a_j,Vmax,Vmin,amax,deltat,N);
    %obj_SCA
    obj_SCA(t) = obj;
    %obj_eUB
    item1 = B*sum(log2(1+gama0./(H^2+vecnorm(q(:,1:N),2).^2)));
    item2 = sum(c1*pow_pos(norms(v(:,1:N),2),3)+c2*inv_pos(norms(v(:,1:N),2)) ... 
    +c2/g^2*quad_over_lin(a,norms(v(:,1:N),2)));
    obj_eUB(t) = item1/item2;
    %obj_ee
    item1 = B*sum(log2(1+gama0./(H^2+vecnorm(q(:,1:N),2).^2)));
    item2 = sum(c1*pow_pos(norms(v(:,1:N),2,1),3)+ c2./norms(v(:,1:N),2).* ...
    (1+(norms(a,2).^2.'-  diag(a.'*v(:,1:N)).^2./norms(v(:,1:N),2).^2.')  /g^2).');
    obj_ee(t) = item1/item2;
    %update
    obj2 = obj;
    epsilonEE = abs(obj2-obj1)/abs(obj2)
    obj1 = obj2;
    q_j = q;
    v_j = v;
    a_j = a;
end
Raver = item1/(N);
Paver = item2/(N);
EE = item1/item2;
%Convergence
figure(1)
title('Convergence of Algorithm 1')
plot(1:t,obj_ee(1:t)/1e3,'-bo','LineWidth',1.3);
hold on,grid on;
plot(1:t,obj_eUB(1:t)/1e3,'--xk','LineWidth',1.3);
plot(1:t,obj_SCA(1:t)/1e3,'--^r','LineWidth',1.3);
xlabel('Iteration')
ylabel('EE(Kbits/Joule)')
legend('1','2','3')
end

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
