function sys=KuramotoLinearSys(A,C,omega,actnode)

  % unweighted graph to weighted graph

n=size(A,1);

%% Compute the desired states

[theta_tgt,omega_tgt]=SolveTgtTheta(A,C,omega,actnode);

%% construct the liear system
A=A.*C;
H=A.*cos(theta_tgt*ones(1,n)-ones(n,1)*theta_tgt');
L2=Adj2Lap(H);
sys.A = -L2;
sys.B = eye(n);
% sys.B=sys.B(:);
sys.theta_tgt=theta_tgt;
sys.omega_tgt=omega_tgt;



end

function L=Adj2Lap(A)

L=-A;

L=L+diag(sum(A,2));

end