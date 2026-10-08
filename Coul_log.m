function lnLambda=Coul_log(chi,tau,a,u,lD)
b0=a.*chi./tau./(3+2*u.^2);
ls2=lD.^2./(1+3./tau./(3+2.*u.^2));
etafit=1+a./sqrt(ls2).*(1+sqrt(1./tau./6));
lnLambda=0.5.*log((b0.^2+etafit.^2.*ls2)./(b0.^2+a.^2));
end
