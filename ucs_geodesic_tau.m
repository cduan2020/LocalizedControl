function [nlist, dlist] = ucs_geodesic_tau(W,node,tau)
% Algorithm 1 (Uniform Cost Search) of the paper
% Calculate the information neighbourhood of information radius tau
% W: upper bound of the information distance matrix
% node: node index of the center of the information neighbourhood
% tau: information radius

if tau==0
    nlist=node;
    dlist=0;
    return;
end

n = size(W,1);
openset = PriorityQueue2(n);

nlist=NaN(1,1000);
dlist=NaN(1,1000);

push(openset, node, 10^-12);

m=1;
gcurrent=0;
while gcurrent<=tau && ~isEmpty(openset)
    
    [u,gu] = pop(openset);
    W(:,u)=0;
    
    [~,neighor,nb_distance]=find(W(u,:));
    
    push(openset, neighor, gu + nb_distance);
        
    
    if size(nlist,2)<m
        nlist=[nlist NaN(1,1000)];
        dlist=[dlist NaN(1,1000)];
    end
    
    nlist(m)=u;
    dlist(m)=gu;   
    gcurrent=gu;
    m=m+1;
end

if gcurrent>tau
    nlist=nlist(1,1:m-2);
    dlist=dlist(1,1:m-2);
else
    nlist=nlist(1,1:m-1);
    dlist=dlist(1,1:m-1);
end

end