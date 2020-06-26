function [theta_tgt]=SolvePostDelta(theta0,V,Y,P,refgen)

% Solve the power flow equation to determine the equilibrium after
% contingencies.

n=length(V);
nonref=setdiff((1:n)',refgen);


fun = @(x) nlsf1(x,V,Y,P,refgen);


max_it=10;
tol=1e-8;
converged=0;
i=0;

x=theta0(nonref);
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    
    %% evaluate Jacobian
    [F,J]=fun(x);

    %% compute update step
    dx = -(J \ F);

    x=x+dx;

    %% check for convergence
    normF = norm(F, inf);
        fprintf('\n%3d        %10.3e', i, normF);
    if normF < tol
        converged = 1;

            fprintf('\nNewton''s method power flow converged in %d iterations.\n', i);

    end
end

if converged==1
    theta_tgt=zeros(n,1);
    theta_tgt(refgen)=0;
    theta_tgt(nonref)=x;
else
    error('Newton Method does not converge !!!');
end


end


function [F,J]=nlsf1(x,Vmag,Y,P,refgen)

n=length(Vmag);
nonref=setdiff((1:n)',refgen);

theta=zeros(n,1);
theta(nonref)=x;

V=Vmag.*exp(1j*theta);

VV=V*V';

I=Y*V;

F=real(V.*conj(I))-P;

F=F(nonref);


J=real(1j*(diag(V.*conj(I))-conj(Y).*VV));

J=J(nonref,nonref);


end
