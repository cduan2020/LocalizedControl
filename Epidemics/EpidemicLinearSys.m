function sys=EpidemicLinearSys(x,Emat,beta,alpha)

% Linearized equation of dynamics

nb=size(Emat,1);

S=x(1:nb);
I=x(nb+1:2*nb);


sys.A = [Emat -diag(beta.*S); 
    diag(beta.*I) Emat-alpha*eye(nb)];

sys.B = -eye(2*nb);


end

