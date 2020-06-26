function Klocal=optcontrol_brain_local(sys,Q,R,L,actnode)

% Local control design for brain network dynamics

ng = size(L,1);
nc = length(R);

Cset0=cell(nc,1);
Cset=cell(nc,1);
Bset=cell(nc,1);
bloc=zeros(nc,1);
Cset_aug0=cell(nc,1);
Cset_aug=cell(nc,1);
AA=cell(nc,1);
BB=cell(nc,1);
QQ=cell(nc,1);
RR=cell(nc,1);
Asize=zeros(nc,1);
Nsize=zeros(nc,1);
ttt=zeros(nc,1);

for i=1:nc
    Cset0{i} = find(L(actnode(i),:));
    
    Cset{i} = find(sum(L((L(actnode(i),:))>0,:),1));
    
    Bset{i}= intersect(Cset{i},actnode);
    bloc(i) = find(Bset{i}==actnode(i));
    
    Cset_aug0{i} = Cset0{i};
    Cset_aug{i} = Cset{i};
    
    AA{i}=sys.A(Cset_aug{i},Cset_aug{i});
    BB{i}=sys.B(Cset_aug{i},Bset{i});
    QQ{i}=Q(Cset_aug{i});
    [~,ia,~]=intersect(actnode,Bset{i}); RR{i}=R(ia);
    Asize(i)=length(Cset0{i});
    Nsize(i)=length(Cset{i});
end


k=cell(nc,1);
parfor i=1:nc
    k{i}=DesignLocalContrl(AA{i},BB{i},QQ{i},RR{i},Cset_aug0{i},Cset_aug{i},bloc(i),ng);
    disp(['system _local_controller_' num2str(i) '_size=' num2str(Asize(i))]);
end
Klocal=vertcat(k{:});



end


