function k=DesignLocalContrl(A,B,Q,R,Cset_aug0,Cset_aug,i,n)

% Local control strategy design by solving the prejected Riccati equation

m=length(R);
[~,~,IB] = intersect(Cset_aug0,Cset_aug);
[Klocal,~,~] = lqr(full(A),full(B),diag(Q),diag(R));
Kex=sparse(m,n);
Kex(:,Cset_aug0)=Klocal(:,IB);
k=Kex(i,:);



end