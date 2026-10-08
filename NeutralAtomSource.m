function NAS=NeutralAtomSource(rd,Rho_d,A,t,Mu)
%A=mI/mi atomic mass = molar wight [g/mol]
% t in seconds =t*rho_s/cs0
NA=6.02e23;
 dt(1)=1;
 
if size(Mu,2)==1 & size(Mu,1)>1 Mu=Mu'; end
dMu=zeros(size(Mu));
for i=1:size(Mu,1)
    for it=2:size(Mu,2)
    dMu(i,it)=Mu(i,it)-Mu(i,it-1);
    dt(i,it)=t(i,it)-t(i,it-1);
    end
end
NAS=-4/3*pi*rd.^3.*Rho_d.*dMu./A.*1000.*NA./dt;