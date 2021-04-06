function dxdt = controlled_brain_eqns(~,x,actnode,Y,alpha,gamma,beta,x_ref,u)

% Brain network dynamics with feedback controllers.

n=size(Y,1);
m=length(actnode);
B=sparse(actnode,1:m,ones(m,1),n,m);


delta = x(1:n);
omega = x(n+1:2*n);
z= x(2*n+1:3*n);

delta_ref = x_ref(1:n);
omega_ref = x_ref(n+1:end);


ddelta_dt =omega;
domega_dt =  -alpha*delta-gamma*delta.^3+beta*Y*delta+B*u;
dz_dt =delta-delta_ref;


dxdt = [ddelta_dt; domega_dt; dz_dt];


end