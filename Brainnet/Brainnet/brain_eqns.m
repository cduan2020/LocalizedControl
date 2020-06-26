function dxdt = brain_eqns(~,x,Y,alpha,gamma,beta)

% Dynamical Euqation of Brian network

n=size(Y,1);

delta = x(1:n);
omega = x(n+1:end);



ddelta_dt =omega;
domega_dt =  -alpha*delta-gamma*delta.^3+beta*Y*delta;

dxdt = [ddelta_dt; domega_dt];


end
