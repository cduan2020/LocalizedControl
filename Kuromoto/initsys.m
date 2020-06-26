function [sys,omega]=initsys(A,C,n,actnode,iter)

% Initialize the control system by calulating the target synchrous state
% and linearized system equation

try
    iter=iter+1;
    omega=(2*pi*rand(n,1)-pi)/5;
    sys=KuramotoLinearSys(A,C,omega,actnode);
    disp(['converged iter=' num2str(iter)]);
catch
    if iter<10
        sys=initsys(A,C,n,actnode,iter);
    else
        error('Still does not converge!!!!!!!!!!!!');
    end
end
end