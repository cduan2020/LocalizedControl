
% Initilized the power grid test case by running the power flow and
% determining the target steady states

%==========================================================================
mpc = ext2int(mpc);
ps = mpc2ps(mpc);
%==========================================================================
[results,success] = runpf(ps); ps = results;

%==========================================================================
% For large-scale system, we may need to get rid of very small generators
% to reduce the stiffness of the differential equation in order to
% numerically simulate the system. 
smallgen=find(ps.gen(:,7)<50);
indices = ps.gen(smallgen,1);
ps.bus(indices,2)=1;
ps.bus(indices,3)=ps.bus(indices,3)-ps.gen(smallgen,2);
ps.bus(indices,4)=ps.bus(indices,4)-ps.gen(smallgen,3);
ps.gen(smallgen,:)=[];
ps.gen_dyn(smallgen,:)=[];
ps.gencost(smallgen,:)=[];
%==========================================================================


%%%%%%%%%%Initialize Machine Variables%%%%%%%%%%%%%%%%%%
swing_index=find(ps.gen(:,1)==find(ps.bus(:,2)==3));  
ref=find(ps.bus(:,2)==3); refgen=find(ps.gen(:,1)==ref);

if isempty(refgen)
    refgen=1;
end

nb=size(ps.bus,1);
ng=size(ps.gen,1);
genbus = ps.gen(:,1);
loadbus = setdiff(1:nb, genbus);

xd=ps.gen_dyn(:,1);
Ra=ps.gen_dyn(:,1)*0;
Yg=1./(Ra+1j*xd);
H=ps.gen_dyn(:,2);
D=ps.gen_dyn(:,3);
wB=2*pi*ps.ref_freq;
M=2*H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


[YredWithGen1,Vref1,delta1,Pm1]=pf2control(ps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


[~,refgen]=max(Pm1);
delta1=delta1-delta1(refgen);

renew_gen=find((0<Pm1)&(Pm1<2));
foss_gen=find(Pm1>=2);


delta3=delta1;
for ratio=1.0:-0.05:0.5
    Pm3=Pm1;
    Pm3(renew_gen)=Pm1(renew_gen)*ratio;
    Pm3(foss_gen)=( sum(Pm1(renew_gen)-Pm3(renew_gen))/sum(Pm3(foss_gen))+1 )*Pm3(foss_gen);
    delta3=SolvePostDelta(delta3,Vref1,YredWithGen1,Pm3,refgen);
end


delta2=delta3;
for ratio=0.5:-0.05:0.2
    Pm2=Pm1;
    Pm2(renew_gen)=Pm1(renew_gen)*ratio;
    Pm2(foss_gen)=( sum(Pm1(renew_gen)-Pm2(renew_gen))/sum(Pm2(foss_gen))+1 )*Pm2(foss_gen);
    delta2=SolvePostDelta(delta2,Vref1,YredWithGen1,Pm2,refgen);
end







