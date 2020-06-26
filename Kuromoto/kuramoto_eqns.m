function dxdt = kuramoto_eqns(~,x,u,A,C,omega)

n=size(A,1);
A=A.*C;
H=A.*sin(x*ones(1,n)-ones(n,1)*x');

dxdt = omega-sum(H,2)+u;


end



