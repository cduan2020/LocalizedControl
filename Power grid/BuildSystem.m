function lure=BuildSystem(YredWithGen,Vref1,wB,M,D,delta1)

% Building the linearized power grid model

n=size(YredWithGen,1);
invM=sparse(diag(1./M));


Pmat=couplePmat(YredWithGen,delta1,Vref1);


lure.A=[sparse(n,n) eye(n); -wB*invM*Pmat -wB*invM*diag(D)];
lure.B=[sparse(n,n); wB*invM];

end