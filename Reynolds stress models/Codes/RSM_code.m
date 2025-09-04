%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REYNOLDS STRESS MODELLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear
format long

load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall friction velocity

viscous=1/395;
nu=viscous;
rho=1;
E=9;
c_mu=0.09;
c1=1.8;
c2=0.6;
c1_p=0.5;
c2_p=0.3;
c1_eps=1.44;
c2_eps=1.92;
sigma_k=1;
sigma_eps=1.3;
kappa=0.41;
urf=0.2;
urf_nu_t=0.5;
S1=10^10;
u_tau=10^(-5);
R=1;
R_max=0.0001;
eps_dns=dns_data(:,2)*1^4/nu;
%Initialize the variables k, epsi, vist, u, u2, v2, w2, uv, f, dudy, duvdy, Uplus


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%equidistant grid
node(1)=0;
node(2)=0.08;
delta=0.01;
n=(1-node(2))/delta+2;
for i=3:n
    node(i)=node(i-1)+delta;
end
face(1)=node(1);
face(2)=node(2)+delta/2;
 for i=3:n
     face(i)=face(i-1)+delta;
 end
 
%To create nodes:the first inner node lies in the log region;y+ = 30 - 100

%To compute delta_yn and delta_ys
delta_yn=delta;
delta_ys=delta;
%intialzation
u=ones(n,1);
uv=ones(n,1);
u2=ones(n,1);
v2=ones(n,1);
w2=ones(n,1);
eps=ones(n,1);
k=ones(n,1);
f=ones(n,1);
dudy=ones(n,1);
duvdy=ones(n,1);
nu_t=ones(n,1);
nu_t_f=ones(n,1);
pk=ones(n,1);
p11=ones(n,1);
p12=ones(n,1);
tau_w=ones(n,1);
s_eps=ones(n,1);
value=ones(n,1);
u_tau=1;
u_tau_old=u_tau;
itr=1;

uexp=interp1(y_dns,u_dns,node,'linear','extrap');
u2exp=interp1(y_dns,u2_dns,node,'linear','extrap');
v2exp=interp1(y_dns,v2_dns,node,'linear','extrap');
w2exp=interp1(y_dns,w2_dns,node,'linear','extrap');
uvexp=interp1(y_dns,uv_dns,node,'linear','extrap');
epsexp=interp1(y_dns,eps_dns,node,'linear','extrap');
u=uexp;
% u2=u2exp;
% v2=v2exp;
% w2=w2exp;
% uv=uvexp;
% eps=epsexp;
u(1)=0;
u2(1)=0;
v2(1)=0;
w2(1)=0;
uv(1)=0;
for i=1:n
        k(i)=(u2(i)+v2(i)+w2(i))/2;
        f(i)=((k(i))^1.5)/(2.55*(node(i))*(eps(i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while R>R_max 
    % compute u_star, tau_w

    itr=itr+1  
    nu_t_old=nu_t;
    u_old=u;
    u2_old=u2;
    v2_old=v2;
    w2_old=w2;
    uv_old=uv;
    eps_old=eps;
    k_old=k;
    u_tau_old= u_tau;

    for i=1:n
        nu_t(i)=rho*c_mu*(k(i)^2)/eps(i);
    end
    
    for i=1:n
        nu_t(i)=0.5*(nu_t(i)+nu_t_old(i));
    end
    
    for i=2:n-1
         nu_t_f(i)=(nu_t(i)+nu_t(i+1))/2;
    end
    
    nu_t_f(n)=nu_t_f(n-1)+(delta/(node(n)-face(n-1)))*(nu_t(n)-nu_t_f(n-1));
    
    % compute the co-effs an,ap,as required for computing 'U'
  
%     for i=2:n-1
%         if i==2
%         an_u(i)=nu/(node(i+1)-node(i));
%         as_u(i)=nu/(node(i)-node(i-1));
%         ap_u(i)=an_u(i)+as_u(i)+u_tau*u_tau/u(i);
%         s(i)=u_tau*u_tau*delta-(1*duvdy(i)*delta);
%         else
%         an_u(i)=nu/(node(i+1)-node(i));
%         as_u(i)=nu/(node(i)-node(i-1));
%         ap_u(i)=an_u(i)+as_u(i);
%         s(i)=1*delta-(1*duvdy(i)*delta);+((u_tau^2/u(i))*delta)
%         end
%      end
 for i=2:n-1
        if i==2
        an_u(i)=nu/(node(i+1)-node(i));
        as_u(i)=nu/(node(i)-node(i-1));
        ap_u(i)=an_u(i)+as_u(i);
        s(i)=1*delta-(1*duvdy(i)*delta);
        else
        an_u(i)=nu/(node(i+1)-node(i));
        as_u(i)=nu/(node(i)-node(i-1));
        ap_u(i)=an_u(i)+as_u(i);
        s(i)=1*delta-(1*duvdy(i)*delta);
        end
 end
   % Wall BC using Wall functionu
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% U EQUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   for i=2:n-1
        u(i)=((an_u(i)*u(i+1))+(as_u(i)*u(i-1))+(s(i)))/ap_u(i); 
   end
   u(n)=u(n-1);
   % relaxation factor
    for i=2:n-1
        u(i)=u_old(i)+urf*(u(i)-u_old(i));
    end

%     u=uexp;
    for i=3:n-1
        
    dudy(1)=(u(2)-u(1))/(node(2)-node(1));
    dudy(2)=u_tau*u_tau/nu;
    dudy(n)=(u(n)-u(n-1))/(node(n)-node(n-1));
    dudy(i)=(u(i+1)-u(i-1))/(node(i+1)-node(i-1));
    
    end
    %   compute u_star, tau_w
    value=(E*u_tau_old*node(2)/nu);
    u_tau=(kappa*u(2))/(log(value));
    u_tau = u_tau_old + 0.5*( u_tau - u_tau_old);
    tau_w=rho*((u_tau)^2);
    
    for i = 1:n
        f (i) = ((k(i))^1.5)/(2.55*(node(i))*(eps(i)));
    end
    
    u2(2)=3.67*(u_tau^2);
    v2(2)=0.83*(u_tau^2);
    w2(2)=2.17*(u_tau^2);
    uv(2)=-1*(u_tau^2);
    eps(2)=u_tau^3/(kappa*node(2));
    % compute eddy viscosity  
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE AND SINK TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        %Production terms
        for i=2:n-1
             p11(i)=-2*rho*uv(i)*dudy(i);
             p12(i)=-rho*v2(i)*dudy(i);
        end
        
           for i=2:n-1
              pk(i)=p11(i)/2;
           end
        
         for i=1:n
             phi_1_w1(i)= c1_p*rho*(eps(i)/k(i))*(v2(i))*(f(i));
             phi_1_w2(i)= c2_p*f(i)*(-c2)*( p11(i)-(2/3)*(p11(i)/2));
             phi_3_w1(i)= phi_1_w1(i);
             phi_12_w2(i)= (-3/2)*c2_p*(f(i))*(-c2*(p12(i)));
         end
        %Pressure strain rate terms
      
      
        %Dissipation terms

        
    %To compute the co-eff an,as required for computing the Reynolds stresses
    % an,as common for u2,v2,w2,uv eqns
      for i=3:n-1
        an(i)=(nu+(nu_t_f(i))/sigma_k)*(1/delta);
        as(i)=(nu+(nu_t_f(i-1))/sigma_k)*(1/delta);
        ap_1(i)=an(i)+as(i)+c1*rho*eps(i)*delta/k(i);
        s_1(i)=p11(i)*delta+((2/3)*c1*rho*eps(i)*delta)-((2/3)*c2*p11(i)*delta)-(2/3)*delta*rho*eps(i)+phi_1_w1(i)*delta+phi_1_w2(i)*delta;
      end
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REYNOLD'S STRESS EQN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %u2
      for i=3:n-1
          u2(i)=(an(i)*(u2(i+1))+as(i)*(u2(i-1))+s_1(i))/ap_1(i);
      end
      u2(n)=u2(n-1);
     for i=3:n-1
         ap_2(i)=an(i)+as(i)+(c1*rho*eps(i)*delta/k(i))+2*c1_p*rho*(eps(i)/k(i))*f(i)*delta;
         s_2(i)=((2/3)*c1*rho*eps(i)*delta)-(2/3)*rho*delta*eps(i)+((2/3)*((c2*pk(i))-(2*c2_p*c2*f(i)*pk(i))))*delta;
         v2(i)=(an(i)*(v2(i+1))+as(i)*(v2(i-1))+s_2(i))/ap_2(i);
     end
     v2(n)=v2(n-1);
     for i=3:n-1
         ap_3(i)=an(i)+as(i)+(c1*rho*eps(i)*delta/k(i));
         s_3(i)=((2/3)*c1*rho*eps(i)*delta)-(2/3)*delta*rho*eps(i)+phi_3_w1(i)*delta+((2/3)*pk(i)*c2*(1 + c2_p*f(i))*delta);
         w2(i)=(an(i)*(w2(i+1))+as(i)*(w2(i-1))+s_3(i))/ap_3(i);
     end
     w2(n)=w2(n-1);
     for i=3:n-1
         ap_12(i)=an(i)+as(i)+(c1*rho*eps(i)*delta/k(i))+c1_p*rho*(eps(i)/k(i))*1.5*f(i)*delta;
         s_12(i)=p12(i)*delta-c2*p12(i)*delta+phi_12_w2(i)*delta;
         uv(i)=(an(i)*uv(i+1)+as(i)*uv(i-1)+s_12(i))/ap_12(i);
     end
      uv(n)=0;
%      for i=3:n
%         u2(i)=u2_old(i) + urf*(u2(i)-u2_old(i));
%         v2(i)=v2_old(i) + urf*(v2(i)-v2_old(i));
%         w2(i)=w2_old(i) + urf*(w2(i)-w2_old(i));
%         uv(i)=uv_old(i) + urf*(uv(i)-uv_old(i));
%      end
    %To compute k, epsilon and f and Uplus
   
    for i=3:n-1
        an_eps(i)=(nu+(nu_t_f(i))/sigma_eps)*(1/delta);
        as_eps(i)=(nu+(nu_t_f(i-1))/sigma_eps)*(1/delta);
        s_eps(i)=c1_eps*rho*eps_old(i)*(p11(i))/(2*k(i));
        ap_eps(i)=an_eps(i)+as_eps(i)+(c2_eps*rho*eps_old(i)*delta/k(i));
        eps(i)=(an_eps(i)*eps(i+1)+as_eps(i)*eps(i-1)+s_eps(i)*delta)/ap_eps(i);
    end
    eps(n)=eps(n-1);
%     for i=3:n-1
%         eps(i)=eps_old(i)+urf*(eps(i)-eps_old(i));
%     end
    for i=1:n
        k(i)=(u2(i)+v2(i)+w2(i))/2;
        f(i)=((k(i))^1.5)/(2.55*(node(i))*(eps(i)));
    end
   for i=2:n-1
        
    duvdy(1)=(uv(2)-uv(1))/(node(2)-node(1));
    duvdy(n)=(uv(n)-uv(n-1))/(node(n)-node(n-1));
    duvdy(i)=(uv(i+1)-uv(i-1))/(node(i+1)-node(i-1));
        
   end
    u=urf*u + ((1-urf)*u_old);
    nu_t=urf_nu_t*nu_t + ((1-urf_nu_t)*nu_t_old);
    k=urf*k + ((1-urf)*k_old);
    u2=urf*u2 + ((1-urf)*u2_old);
    v2=urf*v2 + ((1-urf)*v2_old);
    w2=urf*w2 + ((1-urf)*w2_old);
    uv=urf*uv + ((1-urf)*uv_old);

    eps=urf*eps + ((1-urf)*eps_old);
    %Explicit BC
    R = max([norm(u_old-u),norm(u2-u2_old),norm(v2-v2_old),norm(w2 -w2_old),norm(uv-uv_old),norm(eps_old-eps)]);
    u_plus(1)=0;
    for i=2:n
        y_plus(i)=u_tau*node(i)/nu;
        u_plus(i)=(1/kappa)*log(y_plus(i))+5.2;
    end
    
end

%tau_w=rho*((u_tau)^2);
figure
plot(u_dns,y_dns,'bo'); hold on
plot(u,node,'*r');hold on;
plot(u_plus,node,'*g');
xlabel('u');ylabel('y');
title('Velocity Profile');
legend('DNS','RSTM','Log law');

figure
 %subplot(2,2,1)
 plot(uv_dns,y_dns,'bo'); hold on
plot(uv,node,'*r');
xlabel('uv');ylabel('y');
title('Shear stress');
legend('DNS','RSTM');

%subplot(2,2,2)
figure
plot(u2_dns,y_dns,'bo'); hold on
plot(u2,node,'*r');
xlabel('u2');ylabel('y');
title('Normal stress');
legend('DNS','RSTM');

figure
plot(v2_dns,y_dns,'bo'); hold on
plot(v2,node,'*r');
xlabel('v2');ylabel('y');
title('Normal stress');
legend('DNS','RSTM');

figure
%subplot(2,2,4)
plot(w2_dns,y_dns,'bo'); hold on
plot(w2,node,'*r');
xlabel('w2');ylabel('y');
title('Normal stress');
legend('DNS','RSTM');

figure
plot(eps_dns,y_dns,'bo'); hold on;
plot(eps,node,'*r');
xlabel('eps');ylabel('y');
title('Dissipation');
legend('DNS','RSTM');

figure
k_dns=0.5*(u2_dns+v2_dns+w2_dns);
plot(k_dns,y_dns,'bo'); hold on
plot(k,node,'*r');
xlabel('k');ylabel('y');
title('Turbulent kinetic energy');
legend('DNS','RSTM');