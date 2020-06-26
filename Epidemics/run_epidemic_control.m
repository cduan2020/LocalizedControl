clear;clc;

% Local control design and dynamical simulation of epidemics.

radius = 10;

load('airp_uniq_new.mat');
load('route_new.mat');
load('city_new.mat');

nb=size(airp_uniq_new,1);
nl=size(route_new,1);

A=zeros(nb,nb);
for i=1:nl
    A(route_new(i,1),route_new(i,2))=A(route_new(i,1),route_new(i,2))+route_new(i,3);
end



hs_out=sum(A,2);
hs_in=sum(A,1)';
hs_ave=(hs_out+hs_in);
A(hs_out==0,:) = A(:,hs_out==0)';

city_ind=zeros(size(airp_uniq_new,1),1);
Nairp=zeros(size(airp_uniq_new,1),1);
for i=1:size(airp_uniq_new,1)
    city_ind(i)=airp_uniq_new{i,8};
    Nairp(i)=airp_uniq_new{i,7};
end

for i=1:size(city_uniq,1)
    hh=find(city_ind==i);
    ratio=hs_ave(hh)/sum(hs_ave(hh));
    Nairp(hh)= Nairp(hh).*ratio*1000;
end

hratio=zeros(size(airp_uniq_new,1),1);
for i=1:size(airp_uniq_new,1)
    airp_uniq_new{i,9}=hs_ave(i);
    airp_uniq_new{i,7}=airp_uniq_new{i,7}*1000;
    hratio(i)=airp_uniq_new{i,7}/airp_uniq_new{i,9};
end

S0=Nairp*0.8;
I0=Nairp*0; I0(891)=Nairp(891)*0.01;
beta =  0.1*3.7./Nairp;
alpha = 0.1; % R0=3.7


%% construct coupling matrix

Emat = -diag(sum(A,2)./Nairp) + A./(Nairp*ones(1,nb));

HH= Emat-alpha*eye(nb);

%% weighting matrices and sizes of SIN

Q=[ones(nb,1)*10^-12; ones(nb,1)];
R=[ones(nb,1)*1; ones(nb,1)*100]*5;
invR=1./R;

actnode = 1:nb;

if radius ~= -1
    
    Wup=InfoDistanceUpperbdd(HH);
    l=cell(nb,1);
    for i=1:nb
        [nlist, ~] = ucs_geodesic_k(Wup,i,radius);
        l{i}=sparse(1,nb); l{i}(nlist)=1;
    end
    L=vertcat(l{:});
    
end

%% Simulation


f=@(t,x,u)epidemic_eqns(t,x,u,Emat,beta,alpha);
dt = 1.e-1; T = 600; Num = ceil(T/dt);

xs = zeros(Num+1,2*nb);
us = zeros(Num+1,2*nb);
x=[S0;I0];

x0=x*0;
J=0;
for nt = 1:Num
    xs(nt,:) = x';
    t = (nt-1)*dt
    
    
    if  norm((x-x0)./[Nairp;Nairp],inf)>0.01
        sys=EpidemicLinearSys(x,Emat,beta,alpha);
        Klocal=optcontrol_epidemics_local(sys,Q,R,L,actnode);
        x0=x;
    end
    
    u=-Klocal*x;
    
    
    us(nt,:) = u';
    
    x = rk4_step_ctrl(f,t,x,u,dt);
    x(1:nb) = min(max(x(1:nb),0),Nairp);
    x(nb+1:2*nb) = min(max(x(nb+1:2*nb),0),Nairp);
    
    dJ = ((x-0)'*(Q.*(x-0)) + (u)'*(R.*u))*dt;
    J=J+dJ;
    
end



%% ODE solvers

    function x = rk4_step_ctrl(f,t,x0,u,dt)
        k1 = dt*feval(f,t,x0,u);
        k2 = dt*feval(f,t+0.5*dt,x0 + k1/2,u);
        k3 = dt*feval(f,t+0.5*dt,x0 + k2/2,u);
        k4 = dt*feval(f,t+dt,x0 + k3,u);
        x = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

