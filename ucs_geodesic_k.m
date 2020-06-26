function [nlist, dlist] = ucs_geodesic_k(W,node,k)
% Algorithm 1 (Uniform Cost Search) of the paper
% Calculate the information neighbourhood of size k
% W: upper bound of the information distance matrix
% node: node index of the center of the information neighbourhood
% k: size of the information neighbourhood

if k==1
    nlist=node;
    dlist=0;
    return;
end

n = size(W,1);
openset = PriorityQueue(n);

nlist=NaN(1,k);
dlist=ones(1,k)*30;

push(openset, node, 10^-12);

m=1;
while m<=k && ~isEmpty(openset)
    
    [u,gu] = pop(openset);
    W(:,u)=0;
    
    [~,neighor,nb_distance]=find(W(u,:));
    
    push(openset, neighor, gu + nb_distance);
    
    
    nlist(m)=u;
    dlist(m)=gu;
    m=m+1;
end


nlist = nlist(1:m-1);
dlist = dlist(1:m-1);


end