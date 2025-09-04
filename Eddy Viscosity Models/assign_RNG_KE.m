
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
c_mu=0.085;
%c1_eps=1.44;
c2_eps=1.68;
sigma_k=1;
sigma_eps=0.719;
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
u_RNG=ones(n+1,1);
u_old=u_RNG;
%u_old=ones(n+1,1);
dudy=ones(n+1,1);
dudy_old=ones(n+1,1);
k_RNG=ones(n+1,1);
k_old=ones(n+1,1);
epsilon_RNG=ones(n+1,1);
epsilon_old=ones(n+1,1);
nu_t_RNG=ones(n+1,1);
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
        nu_t_RNG(i)=c_mu*k_RNG(i)^2/epsilon_RNG(i);
    end
    for i=2:n
    fx_n(i)=0.5*(face(i)-face(i-1))/(node(i+1)-node(i));
    fx_s(i)=0.5*(face(i)-face(i-1))/(node(i)-node(i-1));
    end
    nu_t_f(1)=nu_t_RNG(1);
    nu_t_f(n)=nu_t_RNG(n+1);
    for i=2:n
        nu_t_f_n(i)=(fx_n(i)*nu_t_RNG(i+1)+(1-fx_n(i))*nu_t_RNG(i));
        nu_t_f_s(i)=(fx_s(i)*nu_t_RNG(i-1)+(1-fx_s(i))*nu_t_RNG(i));
    end
    
    %Compute U    %----------------------------- NUMERICAL TIP --------------------------
    %Often it can be tricky to start the simulations. They often diverge.
    %It can then be useful to compute the turbulence viscosity from the 
    %mixing-length model for the first few iterations (say 100?).
    %In this way the k & epsilon (or Part II model) equations are de-coupled from U
    %for the first few iterations.
    
    u_RNG(1)=0;
    for i=2:n
         a_s(i)=(nu+nu_t_f_s(i))/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i))/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i);
         u_RNG(i)=(a_s(i)*u_RNG(i-1)+a_n(i)*u_RNG(i+1)+(face(i)-face(i-1)))/a_p(i);
    end
    % u(n)=u(n-1);
     u_RNG(n+1)=u_RNG(n);
     for i=1:n+1
         u_RNG(i)=(u_RNG(i)+u_old(i))*0.5;
     end
        
    %Compute dUdy
    dudy(1)=(u_RNG(2)-u_RNG(1))/(node(2)-node(1));
    dudy(n+1)=(u_RNG(n+1)-u_RNG(n))/(node(n+1)-node(n));
    for i=2:n
        dudy(i)=(u_RNG(i+1)-u_RNG(i-1))/(node(i+1)-node(i-1));
    end
    dudy_f(1)=dudy(1);
    dudy_f(n)=dudy(n+1);
    for i=2:n-1
        dudy_f(i)=(dudy(i+1)-dudy(i))/2;
    end
    %Compute Pk
    for i=1:n+1
        pk_RNG(i,1)=nu_t_RNG(i,1)*(dudy(i,1))^2;
    end
    %Compute k
    k_RNG(1)=0;
   
     for i=2:n
         a_s(i)=(nu+nu_t_f_s(i)/sigma_k)/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i)/sigma_k)/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i);
         b(i)=(pk_RNG(i)-epsilon_RNG(i))*(face(i)-face(i-1));
         k_RNG(i)=(a_s(i)*k_RNG(i-1)+a_n(i)*k_RNG(i+1)+b(i))/a_p(i);
     end
      k_RNG(n+1)=k_RNG(n);
      for i=1:n+1
         k_RNG(i)=(k_RNG(i)+k_old(i))*0.5;
     end

for i=2:n
neta(i)=dudy(i)*k_RNG(i)/epsilon_RNG(i);
end
for i=2:n
    p(i)=neta(i)*(1-(neta(i)/4.38));
    r(i)=1+(0.012*(neta(i))^3);
end
for i=2:n
    c1_eps(i)=1.42-(p(i)/r(i));
end
     for j=2:8
            epsilon_RNG(j)=(nu*2*k_RNG(j))/(node(j)^2);
         end
    for i=9:n
         
         a_s(i)=(nu+nu_t_f_s(i)/sigma_eps)/(node(i)-node(i-1));
         a_n(i)=(nu+nu_t_f_n(i)/sigma_eps)/(node(i+1)-node(i));
         a_p(i)=a_n(i)+a_s(i)+(c2_eps*epsilon_old(i)*(face(i)-face(i-1))/k_RNG(i));
         b(i)=(c1_eps(i)*pk_RNG(i)*epsilon_old(i))*(face(i)-face(i-1))/k_RNG(i);
         epsilon_RNG(i)=(a_s(i)*epsilon_RNG(i-1)+a_n(i)*epsilon_RNG(i+1)+b(i))/a_p(i);
     end
    
    epsilon_RNG(1)=epsilon_RNG(2);
    epsilon_RNG(n+1)=epsilon_RNG(n);
    for i=1:n+1
         epsilon_RNG(i)=(epsilon_RNG(i)+epsilon_old(i))*0.5;
    end
     for i=1:n+1
         shear_stress_RNG(i)= pk_RNG(i)/dudy(i);
     end  
%    
    %Compute residue
    for i=1:n+1
        residue_1=residue_1+abs(u_RNG(i)-u_old(i));
    end
    for i=1:n+1
        residue_2=residue_2+abs(k_RNG(i)-k_old(i));
    end
    for i=1:n+1
        residue_3=residue_3+abs(epsilon_RNG(i)-epsilon_old(i));
    end
    residue_4=max(residue_1,residue_2);
    residue=max(residue_4,residue_3);
    
    epsilon_old=epsilon_RNG;
    k_old=k_RNG;
    u_old=u_RNG;
     
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
    d2k_dy2(i)=(k_RNG(i+1)+k_RNG(i-1)-2*k_RNG(i))/((node(i+1)-node(i-1))^2);
end
d2k_dy2(n+1)=d2k_dy2(n);
dk_dy(1)=(k_RNG(2)-k_RNG(1))/(node(2)-node(1));
dk_dy(2)=(k_RNG(3)-k_RNG(2))/(node(3)-node(2));
d2k_dy2(1)=(dk_dy(2)-dk_dy(1))/(node(2)-node(1));
for i=1:n+1
    k_d_v_RNG(i)=nu*(d2k_dy2(i));
end
for i=1:n+1
    k_d_t_RNG(i)=(nu_t_RNG(i)/sigma_k)*(d2k_dy2(i));
end
for i=1:n+1
    k_d_RNG(i)= k_d_v_RNG(i)+ k_d_t_RNG(i);
end
figure(1)
plot(u_dns,y_dns,'bo');
xlabel('U'); ylabel('y/h'); title('U-velocity');
legend('DNS','Location','Best'); legend boxoff;

figure(2)
plot(y_dns,k_dns,'bo');
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
plot(y_dns,k_d_t_RNG,'r');
xlabel('y/h'); title('Turbulent diffusion of k');
legend('DNS'); legend boxoff

figure(7)
plot(y_dns,dns_data(:,6)/nu,'bo');
hold on
plot(y_dns,k_d_v_RNG,'r');
xlabel('y/h'); title('Viscous diffusion of k');
legend('DNS'); legend boxoff