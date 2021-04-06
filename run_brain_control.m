
% Local control design and dynamical simulation of brain network

clear;clc;
addpath('.\Brainnet');



%% Data Preparation

data=load('Coactivation_matrix');
alpha_p=1935;
alpha_h=2852;
beta=150;
gamma=100;



Y=data.Coactivation_matrix;
n=size(Y,1);


actnode=1:n;

Q=[100*ones(n,1);1e-8*ones(n,1);ones(n,1)];
R=ones(length(actnode),1)*1e-6;


dataMat=beta*Y+alpha_p*eye(n);
Wup=InfoDistanceUpperbdd(dataMat);


%% Build Linear System and SIN

radius=10;
l=cell(n,1);
for i=1:n
    [nlist, ~] = ucs_geodesic_k(Wup,i,radius);
    l{i}=sparse(1,n); l{i}(nlist)=1;
end
L=vertcat(l{:});

x0_h=[0.1*ones(n,1); zeros(n,1);];
x0_p=[0.2*ones(n,1); zeros(n,1); zeros(n,1)];


dt = 1.e-3; T = 5; N = ceil(T/dt);
xs_h = zeros(N+1,2*n);
xs_p = zeros(N+1,3*n);
xs_pc = zeros(N+1,3*n);
f_h = @(t,x)brain_eqns(t,x,Y,alpha_h,gamma,beta);
f_p = @(t,x,x_ref,u)controlled_brain_eqns(t,x,actnode,Y,alpha_p,gamma,beta,x_ref,u);

u0=zeros(length(actnode),1);
x_h=x0_h;
x_p=x0_p;
x_pc=x0_p;

for nt = 1:N
    xs_h(nt,:) = x_h';
    xs_p(nt,:) = x_p';
    xs_pc(nt,:) = x_pc';
    
    t = (nt-1)*dt
    
    sys=BrainLinearSys(x_p,Y,alpha_p,gamma,beta);
    
    
    %     [Klocal1,S,e] = lqr(sys.A,sys.B(:,actnode),diag(Q),diag(R));
    %     K1=full(Klocal1);  % Global control strategy
    
    Klocal=optcontrol_brain_local(sys,Q,R,L,actnode);  % Local control strategy
    
    u = -Klocal*[x_pc(1:n)-x_h(1:n);x_pc(n+1:2*n);x_pc(2*n+1:3*n)];

    x_p = rk4_step_ctrl(f_p,t,x_p,x_h,u0,dt);
    x_pc = rk4_step_ctrl(f_p,t,x_pc,x_h,u,dt);
    x_h = rk4_step_ref(f_h,t,x_h,dt);
end



plot_num=7;

figure(2);
ts = (0:N).'*dt;
subplot(2,1,1);
plot(ts,xs_h(:,plot_num)); title('\delta (radians)'); hold on;
subplot(2,1,2);
plot(ts,xs_h(:,n+plot_num)/(2*pi)); title('\omega (Hz)'); hold on;


figure(2);
ts = (0:N).'*dt;
subplot(2,1,1);
plot(ts,xs_pc(:,plot_num)); title('\delta (radians)'); hold on;
subplot(2,1,2);
plot(ts,xs_pc(:,n+plot_num)/(2*pi)); title('\omega (Hz)'); hold on;


figure(2);
ts = (0:N).'*dt;
subplot(2,1,1);
plot(ts,xs_p(:,plot_num)); title('\delta (radians)'); hold on;
subplot(2,1,2);
plot(ts,xs_p(:,n+plot_num)/(2*pi)); title('\omega (Hz)'); hold on;


%% ODE solvers

function x = rk4_step_ref(f,t,x0,dt)
k1 = dt*feval(f,t,x0);
k2 = dt*feval(f,t+0.5*dt,x0 + k1/2);
k3 = dt*feval(f,t+0.5*dt,x0 + k2/2);
k4 = dt*feval(f,t+dt,x0 + k3);
x = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
end

function x = rk4_step_ctrl(f,t,x0,x_ref,K,dt)
k1 = dt*feval(f,t,x0,x_ref,K);
k2 = dt*feval(f,t+0.5*dt,x0 + k1/2,x_ref,K);
k3 = dt*feval(f,t+0.5*dt,x0 + k2/2,x_ref,K);
k4 = dt*feval(f,t+dt,x0 + k3,x_ref,K);
x = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
end