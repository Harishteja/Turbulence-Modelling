
% DNS data at Re_delta=7890, Re_tau=395 (Moin, Kim & Mansour, PoF, 1999).
% All quantites are normalized by u_tau and nu unless stated otherwise.
% Delta denotes the channel half-width.
clc; close all; clear all

% Read DNS data [half-channel is given (till centerline)]
load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu=1/395;
ustar=1;
rho = 1;
kappa=0.41;
% k-epsilon model constants
c_mu=0.09;
c1_w=(5/9);
c2_w=(3/40);
sigma_k=2;
sigma_w=2;
beta=0.09;
% Part II model constants

% Residual error limit
residue_limit = 10^(-4);
% Under relaxation factor

urf = 0.5;

%Grid (based on DNS data)
%node=y_dns
n=96;
for i=1:n+1
    node(i,1)=y_dns(i);
end
%calculate face values also
face(1,1)= node(1,1);
face(n,1)= node(n+1,1);
% for i=2:n-1
%     face(i,1)=(node(i,1)+node(i+1,1))/2;
% end
 for i=2:n-1
     face(i)=face(i-1)+2*(node(i)-face(i-1));
 end
%Initial conditions for old & new variables (U, dUdy, k, eps, nu_t, residue, ...)
%u=ones(n+1,1);
u_w=ones(n+1,1);
u_old=u_w;
%u_old=ones(n+1,1);
dudy=ones(n+1,1);
dudy_old=ones(n+1,1);
k_w=ones(n+1,1);
k_old=ones(n+1,1);
w=ones(n+1,1);
w_old=ones(n+1,1);
nu_t_w=ones(n+1,1);
nu_t_f=ones(n,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate
residue=1;
z=0;
while residue>residue_limit
    residue_1=0;
    residue_2=0;
    residue_3=0;
    residue_4=0;
    residue=0;
    z=z+1
    %Compute eddy viscosity
    for i=1:n+1
        nu_t_w(i)=k_w(i)/w(i);
    end
    for i=2:n
    fx_n(i)=0.5*(face(i)-face(i-1))/(node(i+1)-node(i));
    fx_s(i)=0.5*(face(i)-face(i-1))/(node(i)-node(i-1));
    end
    nu_t_f(1)=nu_t_w(1);
    nu_t_f(n)=nu_t_w(n+1);
    for i=2:n
        nu_t_f_n(i)=(fx_n(i)*nu_t_w(i+1)+(1-fx_n(i))*nu_t_w(i));
        nu_t_f_s(i)=(fx_s(i)*nu_t_w(i-1)+(1-fx_s(i))*nu_t_w(i));
    end
    
    %Compute U    %----------------------------- NUMERICAL TIP --------------------------
    %Often it can be tricky to start the simulations. They often diverge.
    %It can then be useful to compute the turbulence viscosity from the 
    %mixing-length model for the first few iterations (say 100?).
    %In this way the k & epsilon (or Part II model) equations are de-coupled from U
    %for the first few iterations.
    
    u_w(1)=0;
    for i=2:n
         a_s(i)=(nu+nu_t_f_s(i))/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i))/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i);
         u_w(i)=(a_s(i)*u_w(i-1)+a_n(i)*u_w(i+1)+(face(i)-face(i-1)))/a_p(i);
    end
    % u(n)=u(n-1);
     u_w(n+1)=u_w(n);
     for i=1:n+1
         u_w(i)=(u_w(i)+u_old(i))*0.5;
     end
        
    %Compute dUdy
    dudy(1)=(u_w(2)-u_w(1))/(node(2)-node(1));
    dudy(n+1)=(u_w(n+1)-u_w(n))/(node(n+1)-node(n));
    for i=2:n
        dudy(i)=(u_w(i+1)-u_w(i-1))/(node(i+1)-node(i-1));
    end
    dudy_f(1)=dudy(1);
    dudy_f(n)=dudy(n+1);
    for i=2:n-1
        dudy_f(i)=(dudy(i+1)-dudy(i))/2;
    end
    %Compute Pk
    for i=1:n+1
        pk_w(i,1)=nu_t_w(i,1)*(dudy(i,1))^2;
    end
    %Compute k
    k_w(1)=0;
   
     for i=2:n
         a_s(i)=(nu+nu_t_f_s(i)/sigma_k)/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i)/sigma_k)/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i);
         b(i)=(pk_w(i))*(face(i)-face(i-1));
         k_w(i)=(a_s(i)*k_w(i-1)+a_n(i)*k_w(i+1)+b(i))/(a_p(i)+(0.09*w(i)*(face(i)-face(i-1))));
     end
      k_w(n+1)=k_w(n);
      for i=1:n+1
         k_w(i)=(k_w(i)+k_old(i))*0.5;
      end

     for j=2:8
            w(j)=(nu*6)/((node(j)^2)*c2_w);
         end
    for i=9:n
         
         a_s(i)=(nu+nu_t_f_s(i)/sigma_w)/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i)/sigma_w)/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i)+(c2_w*w_old(i)*(face(i)-face(i-1)));
         b(i)=(c1_w*pk_w(i)*w_old(i))*(face(i)-face(i-1))/k_w(i);
         w(i)=(a_s(i)*w(i-1)+a_n(i)*w(i+1)+b(i))/a_p(i);
     end
    
    w(1)=w(2);
    w(n+1)=w(n);
    for i=1:n+1
         w(i)=(w(i)+w_old(i))*0.5;
    end
    for i=1:n+1
        eps(i)=0.09*k_w(i)*w(i);
    end
     for i=1:n+1
         shear_stress_w(i)= pk_w(i)/dudy(i);
     end
%    
    %Compute residue
    for i=1:n+1
        residue_1=residue_1+abs(u_w(i)-u_old(i));
    end
    for i=1:n+1
        residue_2=residue_2+abs(k_w(i)-k_old(i));
    end
    for i=1:n+1
        residue_3=residue_3+abs(w(i)-w_old(i));
    end
    residue_4=max(residue_1,residue_2);
    residue=max(residue_4,residue_3);
    
    w_old=w;
    k_old=k_w;
    u_old=u_w;
     
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 6 columns in dns_data.dat as below:
%
%      y+         Diss        prod     vel_p_grad   Turb_diff   Visc_diff
%
% Please note that all terms are normalized by ustar^4/nu

k_dns=0.5*(u2_dns+v2_dns+w2_dns);
eps_dns=dns_data(:,2)*ustar^4/nu; % eps is normalized by ustar^4/nu

for i=1:n+1
    nu_t_dns(i)=c_mu*k_dns(i)^2/eps_dns(i);
end
for i=2:n
    d2k_dy2(i)=(k_w(i+1)+k_w(i-1)-2*k_w(i))/((node(i+1)-node(i-1))^2);
end
d2k_dy2(n+1)=d2k_dy2(n);
dk_dy(1)=(k_w(2)-k_w(1))/(node(2)-node(1));
dk_dy(2)=(k_w(3)-k_w(2))/(node(3)-node(2));
d2k_dy2(1)=(dk_dy(2)-dk_dy(1))/(node(2)-node(1));
for i=1:n+1
    k_d_v_w(i)=nu*(d2k_dy2(i));
end
for i=1:n+1
    k_d_t_w(i)=(nu_t_w(i)/sigma_k)*(d2k_dy2(i));
end
for i=1:n+1
    k_d_w(i)= k_d_v_w(i)+ k_d_t_w(i);
end
figure(1)
plot(u_dns,y_dns,'bo');
xlabel('U'); ylabel('y/h'); title('U-velocity');
legend('DNS','Location','Best'); legend boxoff;

figure(2)
plot(y_dns,k_dns,'bo');
hold on
plot(y_dns,k_w,'*')
xlabel('y/h'); ylabel('k'); title('Turbulence kinetic energy');
legend('DNS'); legend boxoff;

figure(3)
plot(y_dns,eps_dns,'bo')
xlabel('y/h'); ylabel('\epsilon'); title('Dissipation rate of k');
legend('DNS'); legend boxoff
    
figure(4)
plot(y_dns,-uv_dns,'bo')
xlabel('y/h'); ylabel('-<uv>'); title('Turbulence shear stress');
legend('DNS'); legend boxoff

figure(5)
plot(y_dns,dns_data(:,3)/nu,'bo');
xlabel('y/h'); ylabel('P_k'); title('Production rate of k');
legend('DNS'); legend boxoff

figure(6)
plot(y_dns,dns_data(:,5)/nu,'bo');
hold on
plot(y_dns,k_d_t_w,'r');
xlabel('y/h'); title('Turbulent diffusion of k');
legend('DNS'); legend boxoff

figure(7)
plot(y_dns,dns_data(:,6)/nu,'bo');
hold on
plot(y_dns,k_d_v_w,'r');
xlabel('y/h'); title('Viscous diffusion of k');
legend('DNS'); legend boxoff