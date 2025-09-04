
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
c1_eps=1.44;
c2_eps=1.92;
sigma_k=1;
sigma_eps=1.3;
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
u=ones(n+1,1);
u_old=u;
%u_old=ones(n+1,1);
dudy=ones(n+1,1);
dudy_old=ones(n+1,1);
k=ones(n+1,1);
k_old=ones(n+1,1);
epsilon=ones(n+1,1);
epsilon_old=ones(n+1,1);
nu_t=ones(n+1,1);
nu_t_f=ones(n,1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate
residue=1;
while residue>residue_limit
    residue_1=0;
    residue_2=0;
    residue_3=0;
    residue_4=0;
    residue=0;
    %Compute eddy viscosity
    for i=1:n+1
        nu_t(i)=c_mu*k(i)^2/epsilon(i);
    end
    for i=2:n
    fx_n(i)=0.5*(face(i)-face(i-1))/(node(i+1)-node(i));
    fx_s(i)=0.5*(face(i)-face(i-1))/(node(i)-node(i-1));
    end
    nu_t_f(1)=nu_t(1);
    nu_t_f(n)=nu_t(n+1);
    for i=2:n
        nu_t_f_n(i)=(fx_n(i)*nu_t(i+1)+(1-fx_n(i))*nu_t(i));
        nu_t_f_s(i)=(fx_s(i)*nu_t(i-1)+(1-fx_s(i))*nu_t(i));
    end
    
    %Compute U    %----------------------------- NUMERICAL TIP --------------------------
    %Often it can be tricky to start the simulations. They often diverge.
    %It can then be useful to compute the turbulence viscosity from the 
    %mixing-length model for the first few iterations (say 100?).
    %In this way the k & epsilon (or Part II model) equations are de-coupled from U
    %for the first few iterations.
    
    u(1)=0;
    for i=2:n
         a_s(i)=(nu+nu_t_f_s(i))/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i))/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i);
         u(i)=(a_s(i)*u(i-1)+a_n(i)*u(i+1)+(face(i)-face(i-1)))/a_p(i);
    end
    % u(n)=u(n-1);
     u(n+1)=u(n);
     for i=1:n+1
         u(i)=(u(i)+u_old(i))*0.5;
     end
        
    %Compute dUdy
    dudy(1)=(u(2)-u(1))/(node(2)-node(1));
    dudy(n+1)=(u(n+1)-u(n))/(node(n+1)-node(n));
    for i=2:n
        dudy(i)=(u(i+1)-u(i-1))/(node(i+1)-node(i-1));
    end
    dudy_f(1)=dudy(1);
    dudy_f(n)=dudy(n+1);
    for i=2:n-1
        dudy_f(i)=(dudy(i+1)-dudy(i))/2;
    end
    %Compute Pk
    for i=1:n+1
        pk(i,1)=nu_t(i,1)*(dudy(i,1))^2;
    end
    %Compute k
    k(1)=0;
   
     for i=2:n
         a_s(i)=(nu+nu_t_f_s(i)/sigma_k)/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i)/sigma_k)/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i);
         b(i)=(pk(i)-epsilon(i))*(face(i)-face(i-1));
         k(i)=(a_s(i)*k(i-1)+a_n(i)*k(i+1)+b(i))/a_p(i);
     end
      k(n+1)=k(n);
      for i=1:n+1
         k(i)=(k(i)+k_old(i))*0.5;
     end

     for j=2:8
            epsilon(j)=(nu*2*k(j))/(node(j)^2);
     end
    for i=9:n
         a_s(i)=(nu+nu_t_f_s(i)/sigma_eps)/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i)/sigma_eps)/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i)+(c2_eps*epsilon_old(i)*(face(i)-face(i-1))/k(i));
         b(i)=(c1_eps*pk(i)*epsilon_old(i))*(face(i)-face(i-1))/k(i);
         epsilon(i)=(a_s(i)*epsilon(i-1)+a_n(i)*epsilon(i+1)+b(i))/a_p(i);
     end
    
    epsilon(1)=epsilon(2);
    epsilon(n+1)=epsilon(n);
    for i=1:n+1
         epsilon(i)=(epsilon(i)+epsilon_old(i))*0.5;
    end
     for i=1:n+1
         shear_stress(i)= pk(i)/dudy(i);
     end   
%    
    %Compute residue
    for i=1:n+1
        residue_1=residue_1+abs(u(i)-u_old(i));
    end
    for i=1:n+1
        residue_2=residue_2+abs(k(i)-k_old(i));
    end
    for i=1:n+1
        residue_3=residue_3+abs(epsilon(i)-epsilon_old(i));
    end
    residue_4=max(residue_1,residue_2);
    residue=max(residue_4,residue_3);
    
    epsilon_old=epsilon;
    k_old=k;
    u_old=u;
     
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
    d2k_dy2(i)=(k(i+1)+k(i-1)-2*k(i))/((node(i+1)-node(i-1))^2);
end
d2k_dy2(n+1)=d2k_dy2(n);
dk_dy(1)=(k(2)-k(1))/(node(2)-node(1));
dk_dy(2)=(k(3)-k(2))/(node(3)-node(2));
d2k_dy2(1)=(dk_dy(2)-dk_dy(1))/(node(2)-node(1));
for i=1:n+1
    k_d_v(i)=nu*(d2k_dy2(i));
end
for i=1:n+1
    k_d_t(i)=(nu_t(i)/sigma_k)*(d2k_dy2(i));
end
for i=1:n+1
    k_d(i)= k_d_v(i)+ k_d_t(i);
end
    k_d_dns=dns_data(:,6)+dns_data(:,5)+dns_data(:,4);
load u_k-w.mat
load u_RNG.mat
load nu_t_w.mat
load nu_t_RNG.mat
load k_RNG.mat
load k-k_w.mat
load epsilon_RNG.mat
load shear_stress_w.mat
load shear_stress_RNG.mat
load pk_RNG.mat
load pk_w.mat
load k_d_v_RNG_f.mat
load k_d_t_RNG_f.mat
load k_d_t_w.mat
load k_d_v_w.mat
load eps.mat
load k_d_w.mat
load k_d_RNG.mat
load nu_t_w.mat
load nu_t_RNG.mat

figure(1)
plot(u_dns,y_dns,'bo');
hold on
plot(u,y_dns,'r','LineWidth',2);
hold on
plot(u_w,y_dns,'g','LineWidth',2);
hold on
plot(u_RNG,y_dns,'b','LineWidth',2);
xlabel('U'); ylabel('y/h'); title('U-velocity');
legend('DNS','k-\epsilon','k-\omega','k-\epsilon (RNG)','Location','Best'); legend boxon;

figure(2)
plot(y_dns,k_dns,'bo');
hold on
plot(y_dns,k,'g','LineWidth',2);
hold on
plot(y_dns,k_RNG,'R','LineWidth',2);
hold on
plot(y_dns,k_w,'*');
xlabel('y/h'); ylabel('k'); title('Turbulence kinetic energy');
legend('DNS','k-\epsilon','k-\epsilon (RNG)','k-\omega'); legend boxon;

figure(3)
plot(y_dns,eps_dns,'bo')
hold on
plot(y_dns,epsilon,'g','LineWidth',2)
hold on
plot(y_dns,epsilon_RNG,'r','LineWidth',2)
hold on
plot(y_dns,eps,'*')
xlabel('y/h'); ylabel('\epsilon'); title('Dissipation rate of k');
legend('DNS','k-\epsilon','k-\epsilon (RNG)','k-\omega'); legend boxon;
    
figure(4)
plot(y_dns,-uv_dns,'bo')
hold on
plot(y_dns,shear_stress,'G','LineWidth',2)
hold on
plot(y_dns,shear_stress_w,'*')
hold on
plot(y_dns,shear_stress_RNG,'r','LineWidth',2)
xlabel('y/h'); ylabel('-<uv>'); title('Turbulence shear stress');
legend('DNS','k-\epsilon','k-\omega','k-\epsilon (RNG)'); legend boxon;

figure(5)
plot(y_dns,dns_data(:,3)/nu,'bo');
hold on
plot(y_dns,pk_RNG,'g','LineWidth',2);
hold on
plot(y_dns,pk_w,'*');
hold on
plot(y_dns,pk,'r','LineWidth',1);
xlabel('y/h'); ylabel('P_k'); title('Production rate of k');
legend('DNS','k-\epsilon(RNG)','k-\omega ','k-\epsilon'); legend boxon;

figure(6)
plot(y_dns,dns_data(:,5)/nu,'bo');
hold on
plot(y_dns,k_d_t,'r','LineWidth',2);
hold on
plot(y_dns,k_d_t_RNG,'g','LineWidth',2);
hold on
plot(y_dns,k_d_t_w,'*');
xlabel('y/h'); title('Turbulent diffusion of k');
legend('DNS','k-\epsilon','k-\epsilon (RNG)','k-\omega'); legend boxon;

figure(7)
plot(y_dns,dns_data(:,6)/nu,'bo');
hold on
plot(y_dns,k_d_v,'r','LineWidth',2);
hold on
plot(y_dns,k_d_v_RNG,'g','LineWidth',2);
hold on
plot(y_dns,k_d_v_w,'*');
xlabel('y/h'); title('Viscous diffusion of k');
legend('DNS','k-\epsilon','k-\epsilon (RNG)','k-\omega'); legend boxon;

figure(8)
plot(y_dns,k_dns,'bo');
hold on
plot(y_dns,-eps_dns,'r','LineWidth',1);
hold on
plot(y_dns,dns_data(:,6)/nu,'g','LineWidth',1);
hold on
plot(y_dns,dns_data(:,5)/nu,'*');
hold on
plot(y_dns,dns_data(:,3)/nu,'--');
xlabel('y/h'); title('TKE-Budget of DNS data');
legend('Turbulence kinetic energy','Dissipation rate of k','Viscous diffusion of k','Turbulent diffusion of k','Production rate of k'); legend boxon

figure(9)
plot(y_dns,k,'bo');
hold on
plot(y_dns,-epsilon,'r','LineWidth',1);
hold on
plot(y_dns,k_d_v,'g','LineWidth',1);
hold on
plot(y_dns,k_d_t,'*');
hold on
plot(y_dns,pk,'--');
xlabel('y/h'); title('TKE-Budget of k-Epsilon Model');
legend('Turbulence kinetic energy','Dissipation rate of k','Viscous diffusion of k','Turbulent diffusion of k','Production rate of k'); legend boxon

figure(10)
plot(y_dns,k_w,'bo');
hold on
plot(y_dns,-eps,'r','LineWidth',1);
hold on
plot(y_dns,k_d_v_w,'g','LineWidth',1);
hold on
plot(y_dns,k_d_t_w,'*');
hold on
plot(y_dns,pk_w,'--');
xlabel('y/h'); title('TKE-Budget of k-omega Model');
legend('Turbulence kinetic energy','Dissipation rate of k','Viscous diffusion of k','Turbulent diffusion of k','Production rate of k'); legend boxon

figure(11)
plot(y_dns,nu_t,'bo');
hold on
plot(y_dns,nu_t_RNG,'r','LineWidth',2);
hold on
plot(y_dns,nu_t_w,'g','LineWidth',2);
xlabel('y/h'); title('Turbulence viscosity');
legend('k-\epsilon','k-\epsilon (RNG)','k-\omega'); legend boxon;