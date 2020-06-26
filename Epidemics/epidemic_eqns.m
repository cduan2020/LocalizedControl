function dxdt = epidemic_eqns(~,x,u,Emat,beta,alpha)

% Dynamical equation of epidemics

nb=size(Emat,1);

S=x(1:nb);
I=x(nb+1:2*nb);

u1=u(1:nb);
u2=u(nb+1:2*nb);

dSdt = -beta.*S.*I+Emat*S-u1;
dIdt = beta.*S.*I-alpha*I+Emat*I-u2;

dxdt = [dSdt ; dIdt];


end