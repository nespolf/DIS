function dydt=eqmotphi(t,y,Ni,M,k,mu,Ti,vphi,vR,R,chi,lnLambda,Te,E,kE,A,Zi,Nn,Tn,Mn,mami)
dydt=NaN.*y;
dydt(1:length(y)/2)=y(length(y)/2+1:end)./R;
%dydt(length(y)/2+1:length(y))= 1./mu.*(zeta_drag(abs(M-y(length(y)/2+1:end)).*sqrt(1./tau),chi,tau,lnLambda).*pi.*mu.^(2/3).*k.*nd.*(M-y(length(y)/2+1:end)).*sqrt(tau)+...
 %   -chi.*Te.*mu.^(1/3).*E.*kE-vphi.*vR./R);
 F=0;
for z=1:length(Zi)
    F=F+zeta_drag(abs(M-y(length(y)/2+1:end))./sqrt(Ti(z)),Zi(z).*chi,Ti(z)/Te,lnLambda(z)).*pi.*mu.^(2/3).*k.*A(z).*Ni(z).*(M-y(length(y)/2+1:end)).*sqrt(Ti(z).*A(1)./A(z));
end
dydt(length(y)/2+1:length(y))= 1./mu.*(F+...
        zeta_drag_n(abs(Mn-y(length(y)/2+1:end))./sqrt(Tn./mami)).*pi.*mu.^(2/3).*k.*mami.*Nn.*(Mn-y(length(y)/2+1:end)).*sqrt(Tn./mami)+...
     -chi.*Te.*mu.^(1/3).*E.*kE)-vphi.*vR./R;
end