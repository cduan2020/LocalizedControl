function Klocal=optcontrol_power_local(lure,Q,R,L,actnode)

% Local control design for power grid dynamics


ng=size(Q,1)/2;
nc=length(actnode);


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

for i=1:nc
    Cset0{i} = find(L(actnode(i),:));
    
    Cset{i} = find(sum(L((L(actnode(i),:))>0,:),1));
    
    
    Bset{i}= intersect(Cset{i},actnode);
    bloc(i) = find(Bset{i}==actnode(i));
    
    Cset_aug0{i} = [Cset0{i} Cset0{i}+ng];
    Cset_aug{i} = [Cset{i} Cset{i}+ng];
    
    AA{i}=lure.A(Cset_aug{i},Cset_aug{i});
    BB{i}=lure.B(Cset_aug{i},Bset{i});
    QQ{i}=Q(Cset_aug{i});
    [~,ia,~]=intersect(actnode,Bset{i}); RR{i}=R(ia);
    Asize(i)=length(Cset0{i});
    Nsize(i)=length(Cset{i});
end


k=cell(nc,1);
for i=1:nc
    k{i}=DesignLocalContrl(AA{i},BB{i},QQ{i},RR{i},Cset_aug0{i},Cset_aug{i},bloc(i),2*ng);
    disp(['system _local_controller_' num2str(i) '_size=' num2str(Asize(i))]);
end
Klocal=vertcat(k{:});


end



