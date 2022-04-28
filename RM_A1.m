function [q,v,a,Raver,Paver,EE] = RM_A1(B,gama0,H,g,c1,c2,q_EE,v_EE,v0,vF,q0,qF,Vmax,Vmin,amax,deltat,T)

%Trajectory optimization to achieve rate-maximization design

%Input
%   B:         Bandwidth
%   gama0:     Reference received signal-to-noise ratio (SNR) at d0 = 1m.
%   H:         Flight Height 
%   g:         Gravitational acceleration
%   c1:        UAV parameter 
%   c2:
%   q_EE:      Initial trajectory
%   v_EE:      Velocity corresponding to q_EE
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
epsilonRM = 1;
th = 1e-3;
%2.SCA and Quadratic Transform
q_j = q_EE;
v_j = v_EE;
t = 0;
obj1 = 0;
while epsilonRM >= th
    %SCA
    t = t+1;
    alpha_j = log2(1+gama0./(H^2+vecnorm(q_j(:,1:N),2).^2)).';
    beta_j = (gama0./(log(2)*(H^2+gama0+vecnorm(q_j(:,1:N),2).^2).*(H^2+vecnorm(q_j(:,1:N),2).^2))).';
 
    cvx_begin quiet
        variable q(2,N+1)
        variable v(2,N+1)
        variable a(2,N)
        expression item1;
        item1 = B*sum(alpha_j-beta_j.* ... 
            (pow_pos(norms(q(:,1:N),2,1),2)-pow_pos(norms(q_j(:,1:N),2,1),2)).');
        maximize item1/1e3  %能解出来
        subject to
            v(:,1) == v0;
            q(:,1) == q0; 
            v(:,N+1) == vF;
            q(:,N+1) == qF;
            norms(v(:,2:N),2) <= Vmax;
            norms(a,2) <= amax;
            pow_pos(norms(v_j(:,1:N),2),2).'+2*diag(v_j(:,1:N).'*(v(:,1:N)-v_j(:,1:N))) >= Vmin^2;
            q(:,2:N+1) == q(:,1:N)+v(:,1:N)*deltat+1/2*a(:,1:N)*deltat^2;
            v(:,2:N+1) == v(:,1:N)+a(:,1:N)*deltat;
    cvx_end
    %update
    obj2 = B*sum(log2(1+gama0./(H^2+vecnorm(q(:,1:N),2).^2)));
    epsilonRM = abs(obj2-obj1)/abs(obj2)
    obj1 = obj2;
    q_j = q;
    v_j = v;
    a_j = a;
end
item1 = B*sum(log2(1+gama0./(H^2+vecnorm(q(:,1:N),2).^2)));
item2 = sum(c1*pow_pos(norms(v(:,1:N),2,1),3)+ c2./norms(v(:,1:N),2).* ...
    (1+(norms(a,2).^2.'-  diag(a.'*v(:,1:N)).^2./norms(v(:,1:N),2).^2.')  /g^2).');
Raver = item1/(N);
Paver = item2/(N);
EE = item1/item2;
end

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
