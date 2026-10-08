function a=NAS_to_EMC3(R,Z,phi,NAS,massrate,A,fnameBsource,Mu)
% R and Z in cm, phi in deg. molI=A_I in g/mol, t in s. NAS in atoms/s
qe=1.6022e-19;
E0=0.356; %eV
NA=6.02e23;
for i=1:size(Mu,1)
    Muend(i)=nanmin(Mu(i,:));
end

if isempty(massrate)
    NAS=NAS./nansum(nansum(NAS))/1000;% normalized to 1 mA total
else
    NAS=NAS./nansum(nansum(NAS)).*massrate/(A.*1e-3)*NA*qe.*(length(Muend)-sum(Muend))./length(Muend); %A*1e-3 molar mass in kg
end
dN=100;
Ntot=round(size(NAS,2)./dN);
NAS1=[]; R1=[]; Z1=[]; phi1=[];
for k=1:size(NAS,1)
for i=1:Ntot
    try
        idx=1+(i-1)*dN:dN.*i;
%     NAS1(i)=sum(NAS(idx));
%     R1(i)=mean(R(idx));
%     Z1(i)=mean(Z(idx));
%     phi1(i)=mean(phi(idx));
    NAS1=[NAS1 sum(NAS(k,idx))];
    R1=[R1 mean(R(k,idx))];
    Z1=[Z1 mean(Z(k,idx))];
    phi1=[phi1 mean(phi(k,idx))];
    catch
    end
end
end
thresh=max(NAS1)./1e4;
R1=R1(NAS1>thresh);
Z1=Z1(NAS1>thresh);
phi1=phi1(NAS1>thresh);
NAS1=NAS1(NAS1>thresh);

a=ones(length(NAS1),8)*0.0;
a(:,1)=E0;
a(:,2)=NAS1;
a(:,3)=R1;
a(:,4)=Z1;
a(:,5)=phi1;
%save(fnameBsource,'a','-ascii')
dlmwrite(fnameBsource,a,'delimiter',' ')