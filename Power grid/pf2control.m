function [YredWithGen,Vref,delta,Pm]=pf2control(ps)

% Convert the target power flow solution to the corresponding control
% strategy of AGC and AVR

nb=size(ps.bus,1);
ng=size(ps.gen,1);
genbus = ps.gen(:,1);
loadbus = setdiff(1:nb, genbus);
ref=find(ps.bus(:,2)==3); refgen=find(ps.gen(:,1)==ref);

if isempty(refgen)
    refgen=1;
end

Pd=ps.bus(:,3)/ps.baseMVA;
Qd=ps.bus(:,4)/ps.baseMVA;
Pg = ps.gen(:,2) / ps.baseMVA; 
Qg = ps.gen(:,3) / ps.baseMVA;
V=ps.bus(:,8);
theta=ps.bus(:,9)*pi/180;
Vsol=V.*exp(1j*theta);


%%%%%%%%%%Initialize Network Variables%%%%%%%%%%%%%%%%%%
MM=full(sparse(genbus,1:ng,ones(ng,1),nb,ng));
E = ((MM'*V + Qg .* xd ./ (MM'*V)) + 1j * (Pg .* xd ./ (MM'*V))) .* exp(1j * MM'*theta);
Vref=abs(E);
delta = angle(E./Vref);
delta=delta-delta(refgen);
Yd=conj((Pd+1j*Qd)./Vsol)./Vsol;
Ywithload=ps.Y+sparse(1:nb,1:nb,Yd); 
YredWithGen= Ywithload(genbus,genbus)-Ywithload(genbus,loadbus)*(Ywithload(loadbus,loadbus)\Ywithload(loadbus,genbus));
Pm=real(diag(Vref.*exp(1j*delta))*conj(YredWithGen)*diag(exp(-1j*delta))*Vref);
norm(Pm-Pg,inf)

end
