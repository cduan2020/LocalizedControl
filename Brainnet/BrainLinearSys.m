function sys=BrainLinearSys(x,Y,alpha,gamma,beta)

% Linearized system of brain network dynamics

n=size(Y,1);
delta=x(1:n);

sys.A=...
    [sparse(n,n), eye(n), sparse(n,n)
    (-alpha)*eye(n)-gamma*diag(delta.^2)+beta*Y,  sparse(n,n),  sparse(n,n)
    eye(n), sparse(n,n), sparse(n,n)];
    
sys.B=...
    [sparse(n,n)
    eye(n)
    sparse(n,n)];

end


