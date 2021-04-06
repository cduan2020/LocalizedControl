function dxdt = swing_eqns(~,x,u,H,wB,M,D,Pm1)

% Dynamical equation of power grid

ng=length(M);

delta = x(1:ng);
omega = x(ng+1:end);

Pe=sum(real((H.*exp(1j*(delta*ones(1,ng)-ones(ng,1)*delta')))),2);
Pm=Pm1+u;
ddelta_dt = omega;
domega_dt =  wB*(Pm - Pe - omega.*D)./M;

dxdt = [ddelta_dt; domega_dt];


end




