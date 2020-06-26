% An example of Algorithm 2 (Gradient-based Greedy Algorithm) of the Paper
% Driverplacment on a BA random network of Laplacian dynamics

clear;
clc;
rng('default');

n=1000;


% Generate network model
A = ba_net('n', n, 'M0', 5, 'M', 5);
C=rand(n,n);
A=A.*C;
A=A-A.*speye(n);
degree = sum(A,2);
L=spdiags(-degree,0,A)-eye(n)*0.1;


% Size of information neighborhood
k=30;

tic;
Wup=InfoDistanceUpperbdd(L);
for i=1:n
    [nlist{i}, dlist{i}] = ucs_geodesic_k(Wup,i,k);
end


%% Driver Placement


ctrlist=nlist;


Wloc=cell(n,1);
for i=1:n
    m=find(ctrlist{i}==i);
    Q=zeros(length(ctrlist{i}),length(ctrlist{i})); Q(m,m)=1;
    Wloc{i}=lyap(L(ctrlist{i},ctrlist{i}),Q);

end



kw = 1; 
ks = 1;


Candidate = 1:n;
Drivenode = [];
Wcest=zeros(n,n);

Eigvec=cell(1,n);
Eigval=zeros(n,1);
for i=1:n
    Eigvec{i}=ones(k,1);
    Eigvec{i}(1)=1;
end
Worstnodes=1:kw;

lamtarget=0.05;
dmax=950;

%% Loop
iter=1;

while (min(Eigval) < lamtarget) && (length(Drivenode)<dmax)
    iter=iter+1;
    
    Bvec = zeros(kw,length(Candidate));
    for j=1:kw
        for i=1:length(Candidate)
            
            if ~isempty(intersect(ctrlist{Candidate(i)},nlist{Worstnodes(j)}(1:k)))
                
                nc=length(ctrlist{Candidate(i)});
                Pc = sparse(1:nc,ctrlist{Candidate(i)},ones(nc,1),nc,n);
                
                Pn = sparse(1:k,nlist{Worstnodes(j)}(1:k),ones(k,1),k,n);
                
                Bvec(j,i) = Eigvec{Worstnodes(j)}'*(Pn*Pc'*Wloc{Candidate(i)}*Pc*Pn')*Eigvec{Worstnodes(j)};
                
            end
        end
    end
    
    [~,ind] = maxk(Bvec,ks);
    chosennode = Candidate(ind);
    
    %==========================================================================================
    
    
    for i=1:length(chosennode)
        Wcest(ctrlist{chosennode(i)},ctrlist{chosennode(i)}) = Wcest(ctrlist{chosennode(i)},ctrlist{chosennode(i)}) + Wloc{chosennode(i)};
    end
    

    Drivenode = union(Drivenode,chosennode);
    Candidate = setdiff(Candidate,chosennode);
    
    
    
    for i=1:n
        [U,D,V]=eig(Wcest(nlist{i}(1:k),nlist{i}(1:k)));
        d=diag(D);
        [mineig,ind] = min(d);
        Eigval(i) = mineig;
        Eigvec{i} = U(:,ind);
    end
    
    
    
    
    [~,Worstnodes] = mink(Eigval,kw);
    disp([num2str(length(Drivenode)), ' have been placed on the network, the estimated controability is ' num2str(Eigval(Worstnodes(1)))]);
    
    
end

% Calculate the true controllability measure
B=eye(n);
B=B(:,Drivenode);
Wc = lyap(L,B*B');
mineig_opt = min(eig(Wc))


% Benchmark against random placement
mineig_vec = zeros(1,1000);

for kk=1:length(mineig_vec)
    kk
    ctrnode = datasample(1:n,dmax,'Replace',false);
    B=eye(n);
    B=B(:,ctrnode);
    Wc = lyap(L,B*B');
    mineig_vec(kk) = min(eig(Wc));
end

