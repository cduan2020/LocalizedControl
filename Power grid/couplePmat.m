function Pmat=couplePmat(YredWithGen,delta,Vref)

% Calcute the coupling matrix P in the linearized power grid dynamics

Y=full(YredWithGen);
V=Vref.*exp(1j*delta);
VV=V*V';
I=Y*V;
Pmat=real(1j*(diag(V.*conj(I))-conj(Y).*VV));


end