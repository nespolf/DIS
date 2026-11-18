function Gamma_vs=vapor_shielding(Te,hsub,mI,kb,qe,ad,epsilon0,me)
%constants
ZW = 15; 
L0 = 1e-31;
RM = 2.5; 
Eev = hsub*mI;
Ted = 1*qe; 
Kew = 1e-13;
IW = 7.86*qe;
ne = 1E19;

%useful quantities
Lambda = log(sqrt(2)*6*pi*(epsilon0*Te)^(3/2)/(ZW^2*sqrt(ne)*qe^2));
kappae = 16*sqrt(2)*epsilon0^2*(Te)^(5/2)/(pi^(3/2)*ZW^2*qe^(3/2)*Lambda*sqrt(me));
alphaL = kappae*L0*mI/(Eev^2);
Gammak = 2*pi*ad^2*kappae*qe*Te/(Eev*sqrt(ad*RM));
Gamma0 = 4*pi*ad*(Ted)^3/(me*Kew*IW*Eev);

%Particule flux
Gamma_vs = Gamma0^(3/5)*Gammak^(2/5)/alphaL^(1/5)*mI;
end
