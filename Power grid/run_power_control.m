% example of power system control 
% with a openly avaiable test system from https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg2000/

clear;clc;
mpc=loadcase('case_ACTIVSg2000');
Init_System;
YredWithGen1=full(YredWithGen1);


delta{1}=delta1; delta{2}=delta2; delta{3}=delta3;
Wup=cell(1,3);
Klocal=cell(1,3);
Asize=cell(1,3);
Nsize=cell(1,3);
tt=cell(1,3);

ng=length(M);
actnode=1:ng;

Q=[ones(ng,1)*10; ones(ng,1)*10];
R=ones(length(actnode),1);
invR=1./R;

radius=10;
sys = cell(3,1);
for i=1:3
    sys{i}=BuildSystem(YredWithGen1,Vref1,wB,M,D,delta{i});
    Wup=InfoDistanceUpperbdd(RecoverCouple_full(sys{i}.A,ng));
    l=cell(ng,1);
    for ii=1:ng
        [nlist, ~] = ucs_geodesic_k(Wup,ii,radius);
        l{ii}=sparse(1,ng); l{ii}(nlist)=1;
    end
    L=vertcat(l{:});
    Klocal{i}=optcontrol_power_local(sys{i},Q,R,L,actnode);
end


my_ode = @rk4_step;
H=conj(YredWithGen1).*(Vref1*Vref1');
f = @(t,x,u,Pm1)swing_eqns(t,x,u,H,wB,M,D,Pm1);
dt = 1.e-3; T = 10; N = ceil(T/dt);

omega1=zeros(ng,1);
x = [delta1; omega1]; xs = zeros(N+1,2*ng);  us = zeros(N+1,ng);
t=0;
J=0;
for n = 1:N
    
    if t<1
        xs(n,:) = x;
        t = (n-1)*dt;
        x_target=[delta1; zeros(ng,1)];
        u=-Klocal{1}*(x-x_target);
        us(n,:) = u+Pm1;
        x = my_ode(f,t,x,u,Pm1,dt);
    elseif (1<=t) && (t<4)
        xs(n,:) = x;
        t = (n-1)*dt;
        x_target=[delta2; zeros(ng,1)];
        u=-Klocal{2}*(x-x_target);
        us(n,:) = u+Pm2;
        x = my_ode(f,t,x,u,Pm2,dt);
    elseif (4<=t) && (t<7)
        xs(n,:) = x;
        t = (n-1)*dt;
        x_target=[delta3; zeros(ng,1)];
        u=-Klocal{3}*(x-x_target);
        us(n,:) = u+Pm3;
        x = my_ode(f,t,x,u,Pm3,dt);
    else
        xs(n,:) = x;
        t = (n-1)*dt;
        x_target=[delta1; zeros(ng,1)];
        u=-Klocal{1}*(x-x_target);
        us(n,:) = u+Pm1;
        x = my_ode(f,t,x,u,Pm1,dt);
    end
    dJ = ((x-x_target)'*(Q.*(x-x_target)) + (u)'*(R.*u))*dt;
    J=J+dJ;
    disp(['radius = ' num2str(radius) ', time = ' num2str(t) ', J = ' num2str(J) ', dJ = ' num2str(dJ)]);
    
end

figure;
ts = (0:N).'*dt;
subplot(2,1,1);
plot(ts,xs(:,1:ng)); title('\delta (radians)');
subplot(2,1,2);
plot(ts,xs(:,ng+1:end)/(2*pi)+60); title('\omega (Hz)');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE solvers
function x = rk4_step(f,t,x0,u,Pm,dt)
  k1 = dt*feval(f,t,x0,u,Pm);
  k2 = dt*feval(f,t+0.5*dt,x0 + k1/2,u,Pm);
  k3 = dt*feval(f,t+0.5*dt,x0 + k2/2,u,Pm);
  k4 = dt*feval(f,t+dt,x0 + k3,u,Pm);
  x = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
end

% Coupling strength matrix
function P=RecoverCouple_full(A,n)

T11=abs(A(1:n,1:n));
T12=abs(A(1:n,n+1:end));
T21=abs(A(n+1:end,1:n));
T22=abs(A(n+1:end,n+1:end));
P=max(max(max(T11,T12),T21),T22);
P=max(abs(P),abs(P)');

end

