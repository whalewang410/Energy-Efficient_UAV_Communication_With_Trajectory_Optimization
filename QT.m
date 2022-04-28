function [q,v,a,obj2] = QT(g,B,alpha_j,beta_j,c1,c2,v0,vF,q0,qF,v_j,q_j,a_j,Vmax,Vmin,amax,deltat,N)

%Quadratic transform method to solve fractional programming

%1.Initialize
item1 = B*sum(alpha_j);
item2 = sum(c1*vecnorm(v_j(:,1:end-1)).^3+c2./vecnorm(v_j(:,1:end-1)) ... 
    +c2*(vecnorm(a_j).^2)./(g^2*vecnorm(v_j(:,1:end-1))));   %1-(N-1)的速度耗能
obj1 = 0;
epsilon_QT = 1;
th = 1e-3;

%AO
while epsilon_QT >= th
    y_t = sqrt(item1)/item2;
    cvx_begin quiet
        variable q(2,N+1)
        variable v(2,N+1)
        variable a(2,N)
        variable tau(N);
        expression item1;
        expression item2;
        item1 = B*sum(alpha_j-beta_j.* ... 
            (pow_pos(norms(q(:,1:N),2,1),2)-pow_pos(norms(q_j(:,1:N),2,1),2)).');
        item2 = sum(c1*pow_pos(norms(v(:,1:N),2,1),3).'+c2*inv_pos(tau) ... 
    +c2/g^2*quad_over_lin(a,tau.').');

        maximize ((2*y_t*sqrt(item1)-y_t^2*item2)*1e5)
        subject to
            v(:,1) == v0;
            q(:,1) == q0; 
            v(:,N+1) == vF;
            q(:,N+1) == qF;
            norms(v(:,2:N),2) <= Vmax;
            norms(a,2) <= amax;
            tau >= Vmin;
            pow_pos(norms(v_j(:,1:N),2),2).'+2*diag(v_j(:,1:N).'*(v(:,1:N)-v_j(:,1:N))) >= tau.^2;
            q(:,2:N+1) == q(:,1:N)+v(:,1:N)*deltat+1/2*a(:,1:N)*deltat^2;
            v(:,2:N+1) == v(:,1:N)+a(:,1:N)*deltat;
    cvx_end
    obj2 = cvx_optval/1e5;
    epsilon_QT = abs(obj2-obj1)/abs(obj2);
    obj1 = obj2;
end
end

                
                
                
                
                
                
                
                
                
                
                
                
                