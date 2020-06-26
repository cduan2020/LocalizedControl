clear;clc;

% Local control design and dynamical simulation of Kuramoto oscillators

%% generate random networks
n=1000;
A = ba_net('n', n, 'M0', 5, 'M', 5);
n=size(A,1);
cvalue=1;
C=cvalue*rand(n,n).*A;


%% generate kuramoto model


actnode=1:n; % list of driver nodes

[sys,omega]=initsys(A,C,n,actnode,0); % calculate the sync state and the linearized system model

Np=1;
xs=2*pi*rand(n,Np);


Q=ones(n,1)*5;
R=ones(length(actnode),1);

f=@(t,x,u)kuramoto_eqns(t,x,u,A,C,omega);
dt = 1.e-1; T = 10; Num = ceil(T/dt);

%% construct the information structure graph

dataMat=sys.A;
Wup=InfoDistanceUpperbdd(dataMat);
radius=10;
l=cell(n,1);
for i=1:n
    [nlist, ~] = ucs_geodesic_k(Wup,i,radius);
    l{i}=sparse(1,n); l{i}(nlist)=1;
end
L=vertcat(l{:});

%% Central Control Design

[K1,S,e] = lqr(sys.A,sys.B(:,actnode),diag(Q),diag(R));


%% Local Control Design

Klocal=optcontrol_brain_local(sys,Q,R,L,actnode);


%% Simulation

Ktest=Klocal;
% Ktest=K1;
x=xs(:,1);
J=0;
xvec2=zeros(n,Num);
for nt = 1:Num
    t = (nt-1)*dt;
    u=-Ktest*sin(x-sys.theta_tgt-sys.omega_tgt*t)+sys.omega_tgt-omega(actnode);
    x =wrapToPi( rk4_step_ctrl(f,t,x,sys.B(:,actnode)*u,dt));
    x_tgt=wrapToPi(sys.theta_tgt+sys.omega_tgt*t);
    dJ=(x-x_tgt)'*(Q.*(x-x_tgt))*dt+(u-(sys.omega_tgt-omega(actnode)))'*(R.*(u-(sys.omega_tgt-omega(actnode))))*dt;
    J=J+dJ;
    xvec2(:,nt)=x;
end




ts = (0:Num-1).'*dt;
plot(ts,xvec2);


function x = rk4_step_ctrl(f,t,x0,u,dt)
k1 = dt*feval(f,t,x0,u);
k2 = dt*feval(f,t+0.5*dt,x0 + k1/2,u);
k3 = dt*feval(f,t+0.5*dt,x0 + k2/2,u);
k4 = dt*feval(f,t+dt,x0 + k3,u);
x = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
end



