function zeta=zeta_drag(x,chi,tau,lnLambda)
% zeta=zeta_drag(x,chi,tau,lnLambda)
% from [Smirnov PPCF 2007]
xxp=x.^2+chi./tau;
xxm=x.^2-chi./tau;
xp=x+sqrt(-chi./tau);
xm=x-sqrt(-chi./tau);
if chi>=0
zeta_coll=1./(2*x.^2).*(1/sqrt(pi).*(1+2*xxp).*exp(-x.^2)+x.*(1+2.*xxp-1./(2.*x.^2).*(1-2.*xxm)).*erf(x) );
else
zeta_coll=1./(4*x.^2).*(1/sqrt(pi).*((1+2.*x.^2+(1-2.*x.^2)./x.*sqrt(-chi./tau)).*exp(-xp.^2)+...
    (1+2.*x.^2-(1-2.*x.^2)./x.*sqrt(-chi./tau)).*exp(-xm.^2))...
    +x.*(1+2.*xxp-1./(2.*x.^2).*(1-2.*xxm)).*(erf(xp)+erf(xm)));   
end
zeta_orb=2.*lnLambda.*(chi./tau).^2.*(erf(x)-2.*x./sqrt(pi).*exp(-x.^2))./(2.*x.^3);
zeta=zeta_coll+zeta_orb;
if isnan(zeta) zeta=0; end
end