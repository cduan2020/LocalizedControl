function [theta_tgt,omega_tgt]=SolveTgtTheta(A,C,omega,actnode)

% Calculate the target synchrous state

A=A.*C;
n=size(A,1);
freenode=setdiff(1:n,actnode);
omega_tgt=mean(omega(freenode));


if isempty(freenode)
    theta_tgt=zeros(n,1);
    omega_tgt=0;
    return;
end


x0=zeros(length(freenode),1);

fun = @(x) nlsf1(x,n,A,omega,actnode);
% options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true,'MaxIter',10000);
% [x,F,exitflag,output,JAC] = fsolve(fun,x0,options);

[x, resnorm, F, exitflag, output, jacob] = newtonraphson(fun, x0);

if exitflag==1
    theta_tgt=zeros(n,1);
    theta_tgt(actnode)=0;
    theta_tgt(freenode)=x;
else
    error('Newton Method does not converge !!!');
end


end




