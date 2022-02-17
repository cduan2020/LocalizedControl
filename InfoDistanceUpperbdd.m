function Wup=InfoDistanceUpperbdd(dataMat)


% Building the the matrix Wup, the upper bound for information distances,
% from the system coupling matrix dataMat

% parameters of the sub-exponential function v()
alpha=1;
beta=0.9999;
s=1.0001;

n=size(dataMat,1);
x=0:0.01:30;
y=exp(alpha*x.^(beta)).*(1+x).^s;

HH=abs(dataMat);


c=full(max(max((HH))));


[i,j,v] = find(HH);
v2d = interp1(y,x,c./v) + 1e-12;

Wup = sparse(i,j,v2d,n,n);
R = spones(Wup);

Wup=(R.*R').*(min(Wup,Wup'))+(1-R.*R').*(max(Wup,Wup'));
Wup(1:n+1:n*n)=0;


end
