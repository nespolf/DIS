function zeta=zeta_drag_n(u)
% zeta=zeta_drag_n(u)
% from [Smirnov PPCF 2007]
zeta=1./u.*(1./sqrt(pi).*(u+1./2./u).*exp(-u.^2)+(1+u.^2-1./4./u.^2).*erf(u));
if isnan(zeta) zeta=0; end
end