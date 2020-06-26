function averge_neighborhood_rate = average_k_rate(L,k)

% Calculate the average Lth neighbor reduction rate

n=size(L,1);
Wup=InfoDistanceUpperbdd(L);
nlist=cell(1,n);
dlist=cell(1,n);
neighbor_rate=zeros(n,k);


alpha=1;
beta=0.9999;
s=1.0001;
y= @(x) exp(alpha*x.^(beta)).*(1+x).^s;

absL = abs(L);
c = full(max(max(absL)));
d = max(max(max(max(absL,[],1)',max(absL,[],2)),1e-12),1);

parfor i=1:n
    [nlist{i}, dlist{i}] = ucs_geodesic_k2(Wup,i,k);
    neighbor_rate(i,:) = c./y(dlist{i})./d(i);
end

averge_neighborhood_rate = mean(neighbor_rate,1);

end