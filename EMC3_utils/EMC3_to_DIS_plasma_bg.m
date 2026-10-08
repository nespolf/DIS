clear all
try
tic;
fnameout=['DIS_pbg_EMC3.mat'];
inputgeofile='/p/emc3-eirene/fnespoli/LHD_grids/360-036/input.geo';
gridfile='/p/emc3-eirene/fnespoli/LHD_grids/360-036/open/lhd_gr.txt';
kmagfile='/p/emc3-eirene/fnespoli/kmag/output.txt';

% inputgeofile='/Users/fnespoli/EMC3_post/LHD_grids/360-036/input.geo';
% gridfile='/Users/fnespoli/EMC3_post/LHD_grids/360-036/open/lhd_gr.txt';
% kmagfile='/p/emc3-eirene/fnespoli/kmag/output.txt';

% %% read grid info from stdout
% disp('read grid info...')
% fname=dir('stdout*');
% fid=fopen(fname(1).name);
% for i=1:15
%     a=fgetl(fid);
% end
% fclose(fid);
% NZ=str2num(a(8:end));
% range=['17:',num2str(17+NZ-1)];
% out=readmatrix(fname(1).name,'FileType','text','range',range);
% out=out(:,1:4);
% IZ=out(:,1)+1;
% Nr=out(:,2);
% Nth=out(:,3);
% Nphi=out(:,4);
% clear out Mesh_p_os Grid_p_os
% Mesh_p_os=0; Grid_p_os=0;
% for iz=1:length(IZ)-1
%     Mesh_p_os(iz+1)=(Nr(iz)-1)*(Nth(iz)-1)*(Nphi(iz)-1)+Mesh_p_os(iz);
%     Grid_p_os(iz+1)=(Nr(iz))*(Nth(iz))*(Nphi(iz))+Grid_p_os(iz);
% end
%% read grid info from input.geo
disp('read grid info...')
fid=fopen(inputgeofile);
for i=1:2
    a=fgetl(fid);
end
NZ=str2num(a); IZ=1:NZ;
clear Nr Nth Nphi
for i=1:NZ
    a=str2num(fgetl(fid));
    Nr(i)=a(1);
    Nth(i)=a(2);
    Nphi(i)=a(3);
end
fclose(fid);
Nr=Nr'; Nth=Nth'; Nphi=Nphi';
clear out Mesh_p_os Grid_p_os
Mesh_p_os=0; Grid_p_os=0;
for iz=1:length(IZ)-1
    Mesh_p_os(iz+1)=(Nr(iz)-1)*(Nth(iz)-1)*(Nphi(iz)-1)+Mesh_p_os(iz);
    Grid_p_os(iz+1)=(Nr(iz))*(Nth(iz))*(Nphi(iz))+Grid_p_os(iz);
end
%% read the grid
disp('read grid ...')

clear R Z PHI
fid=fopen(gridfile);
[a]=fgetl(fid);
n=str2num(a);
i=1; R1=[]; Z1=[];  step=1; iz=1; Nzones=11;
while (~feof(fid))
    
    a=str2num(fgetl(fid));
    if length(a)==1 && step==1
        Phi(i,iz)=a; step=2;
    elseif step==2 && length(a)==6 && mean(a)<200
        Z1=[Z1 a]; step=3;
    elseif step==2 && length(a)==6
        R1=[R1 a];
    elseif (step==2 && length(a)<6)
        R1=[R1 a]; step=3;
    elseif step==3 && length(a)==6
        Z1=[Z1 a];
    elseif step==3 && length(a)<6
        Z1=[Z1 a];
        R(iz).d(i,:)=R1;
        Z(iz).d(i,:)=Z1;
        PHI(iz).d(i,:)=ones(size(R1)).*Phi(i,iz);
        i=i+1;
        R1=[]; Z1=[];
        step=1;
        if i>n(3) i=1; iz=iz+1; end
    end
    
end
%%
RG=[]; ZG=[]; PHIG=[];

Mesh_p_os(1)=0;
for iz=1:length(R)
    a=R(iz).d';
    RG=[RG a(:)'];
    b=Z(iz).d';
    ZG=[ZG b(:)'];
    c=PHI(iz).d';
    PHIG=[PHIG c(:)'];
end


%% read cell geo
disp('read cell geo...')

clear IDCELL
out=readmatrix('CELL_GEO','FileType','text','NumHeaderLines',1)';
IDCELL=out(:);
fid=fopen('CELL_GEO');
Ncell=str2num(fgetl(fid));
% if isnan(IDCELL)%case from stellar with CELL_GEO on a single line
% IDCELL=str2num(fgetl(fid))+1;
% end
fclose(fid);
if isnan(IDCELL)
fid=fopen('CELL_GEO');
formatSpec = '%f';
A = fscanf(fid,formatSpec);
fclose(fid)
Ncell=A(1:3);
IDCELL=A(4:end);
clear A
end

%% read Temperature
disp('read Temperature...')

clear out
out=readmatrix('TE_TI','FileType','text')';
out=out(:);
idx=find(isnan(out));
out=out(isnan(out)<1);
TE=out(1:idx(1)-1);
TI=out(idx(1):end);
%% read density 
disp('read density...')
try
clear out N idx1 idx2
out=readmatrix('DENSITY','FileType','text')';

out=out(:);
idx=find(isnan(out));
idxprev=0; k=0; k1=0; idx1(1)=0;
for i=1:length(idx)
     if idx(i)-idxprev>Ncell(2)
         k=k+1;
         idx1(k)=idx(i); 
         idxprev=idx(i); 
     end
     if i>1 && idx(i)-idx(i-1)>Ncell(2)
         k1=k1+1;
         idx2(k)=idx(i-1); 
     end
end

for i=1:length(idx1)
N(i,:)=out(idx2(i)+1:idx1(i)-1);
end
catch
 
    %% read density 2nd method
disp('read Density 2nd method...')
clear out N A
out=readmatrix('DENSITY','FileType','text','NumHeaderLines',0)';
linnum=1:size(out,2);
linnum=linnum(isnan(sum(out,1)));
if linnum(end)<size(out,2)
    linnum(end+1)=size(out,2)+1;
end
for i=1:length(linnum)-1
    A=out(:,linnum(i)+1:linnum(i+1)-1);
    N(i,:)=A(:);
end
if size(N,2)~=Ncell(2) 
     disp('read Density 3nd method...')
clear out N 
out=readmatrix('DENSITY','FileType','text','NumHeaderLines',0)';
out=out(:);
out=out(isnan(out)<1); 
Ztot=floor(length(out)/Ncell(2));
for i=1:Ztot
N(Ztot,:)=out((i-1).*Ncell(2)+1+i:i.*Ncell(2)+i);
end
end

end
%% read neutral density
disp('read neutral density...')

clear out
out=readmatrix('DENSITY_A','FileType','text')';
out=out(:);
out=out(isnan(out)<1);
Nn=out(:);
%% read neutral Temperature
disp('read neutral Temperature...')

clear out
out=readmatrix('TEMPERATURE_A','FileType','text')';
out=out(:);
out=out(isnan(out)<1);
Tn=out(:);


%% read mach number
disp('read mach number...')

clear out M
out=readmatrix('MACH_NUMBER','FileType','text')';
out=out(isnan(out)<1);
M=out(:);

%% read impurity velocity
disp('read impurity velocity...')

if size(N,1)>1
    clear out MI idx idx1 idx2
    out=readmatrix('IMPURITY_VELOCITY_2','FileType','text')';
    
    out=out(:);
    idx=find(isnan(out));
    idxprev=0; k=0; k1=0; idx1(1)=0;
    for i=1:length(idx)
        if idx(i)-idxprev>Ncell(2)
            k=k+1;
            idx1(k)=idx(i);
            idxprev=idx(i);
        end
        if i>1 && idx(i)-idx(i-1)>Ncell(2)
            k1=k1+1;
            idx2(k)=idx(i-1);
        end
    end
    
    for i=1:length(idx1)
        MI(i,:)=out(idx2(i)+1:idx1(i)-1);
    end
end

%% read impurity radiation
try
disp('read impurity radiation...')

clear out IR
out=readmatrix('IMP_RADIATION','FileType','text')';
IR=out(:);
IR=IR(~isnan(IR));
catch
    IR=0;
end
%[W/cc]
 %% read plasma radiation
 try
disp('read plasma radiation...')

clear out PR
out=readmatrix('RADIATION_1','FileType','text')';
PR=out(:);
PR=PR(~isnan(PR));
 catch
     PR=0
 end
%% total radiation
RAD=(PR-IR).*1e6; %[W/m3]
clear PR IR
%%
disp('remapping the cells...')

clear IG Np Rp Zp PHIp IG1 I1 I2 I3 I4 I5 I6 I7 I8

for iz=1:length(R)
    for K=0:Nphi(iz)-2
        lll=1;
        for J=0:Nth(iz)-2
            for I=0:Nr(iz)-2
                IG1(iz).d(K+1,lll)=I+J*(Nr(iz)-1)+K*(Nr(iz)-1)*(Nth(iz)-1)+Mesh_p_os(iz);
                I1(iz).d(K+1,lll)= I  + J  *(Nr(iz))+K*(Nr(iz))*(Nth(iz))+Grid_p_os(iz);
                I2(iz).d(K+1,lll)= I+1+ J  *Nr(iz) + K*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                I3(iz).d(K+1,lll)= I+1+ (J+1)*Nr(iz) + K*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                I4(iz).d(K+1,lll)= I  + (J+1)*Nr(iz) + K*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                I5(iz).d(K+1,lll)= I  + J  *Nr(iz) + (K+1)*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                I6(iz).d(K+1,lll)= I+1+ J  *Nr(iz) + (K+1)*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                I7(iz).d(K+1,lll)= I+1+ (J+1)*Nr(iz) + (K+1)*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                I8(iz).d(K+1,lll)= I  + (J+1)*Nr(iz) + (K+1)*Nr(iz)*Nth(iz) +Grid_p_os(iz);
                lll=lll+1;
            end
            
        end
        
    end
end

IG=[]; I11=[]; I21=[]; I31=[]; I41=[]; I51=[]; I61=[]; I71=[]; I81=[];
for iz=1:length(IG1)
    a=IG1(iz).d'; IG=[IG a(:)'];
    a=I1(iz).d'; I11=[I11 a(:)'];
    a=I2(iz).d'; I21=[I21 a(:)'];
    a=I3(iz).d'; I31=[I31 a(:)'];
    a=I4(iz).d'; I41=[I41 a(:)'];
    a=I5(iz).d'; I51=[I51 a(:)'];
    a=I6(iz).d'; I61=[I61 a(:)'];
    a=I7(iz).d'; I71=[I71 a(:)'];
    a=I8(iz).d'; I81=[I81 a(:)'];
end
IG=IG+1;
I11=I11+1;
I21=I21+1;
I31=I31+1;
I41=I41+1;
I51=I51+1;
I61=I61+1;
I71=I71+1;
I81=I81+1;

%%
clear R1 Z1 PHI1 Np Rp Zp Phip Mp Tep Tip  Npn Rpn Zpn PHIpn Tpn mask maskn RADp
Ztot=size(N,1);
k=1; k1=1;
for i=1:length(IG)%Ncell(2)%length(IC)
    IC=int32(IDCELL(IG(i)));
    idx=[I11(i) I21(i) I31(i) I41(i) I51(i) I61(i) I71(i) I81(i)];
    R1(i)=mean(RG(idx));
    Z1(i)=mean(ZG(idx));
    PHI1(i)=mean(PHIG(idx));
    if IC<=Ncell(2)
        Np(:,k)=N(:,IC);
        if Ztot>1
        MIp(:,k)=MI(:,IC);
        end
        Mp(k)=M(IC);
        Tep(k)=TE(IC);
        Tip(k)=TI(IC);
        RADp(k)=RAD(IC);
        Rp(k)=R1(i);  Zp(k)=Z1(i); PHIp(k)=PHI1(i); mask(i)=1;
        k=k+1;
    else
        R1(i)=mean(RG(idx));
        Z1(i)=mean(ZG(idx)); PHI1(i)=mean(PHIG(idx));  mask(i)=NaN;
    end
    if IC<=Ncell(3)
        Npn(k1)=Nn(IC); Tpn(k1)=Tn(IC);
        Rpn(k1)=R1(i);  Zpn(k1)=Z1(i); PHIpn(k1)=PHI1(i);
        k1=k1+1;
    else
    end
end


%%
clear R Z RG ZG PHIG N TE TI M MI I1 I2 I3 I4 I5 I6 I7 I8  I11 I21 I31 I41 I51 I61 I71 I81 IG IDCELL Nn Tn
Rp=Rp/100; Zp=Zp/100; R1=R1/100; Z1=Z1/100; Np=Np.*1e6;
Rpn=Rpn/100; Zpn=Zpn/100; Npn=Npn.*1e6;


Nr=498; Nth=502;
%Nr=250; Nth=250;

R=linspace(min(Rp),max(Rp),Nr);
Z=linspace(min(Zp),max(Zp),Nth);
phi=linspace(min(PHIp),max(PHIp),64);

[R Z]=meshgrid(R,Z);

clear fieldN fieldTe fieldTi fieldM fieldNn fieldTn fieldMI fieldPrad
A=Rpn(PHIpn==min(PHIpn)); B=Zpn(PHIpn==min(PHIpn));
k=boundary(A(:),B(:));
A=A(k); B=B(k);
disp('interpolating mask')

maskn=double(inpolygon(R,Z,A,B));
maskn(maskn==0)=NaN;
%%
disp('### Mask')
disp('creating scatteredinterpolant...')
Mask=scatteredInterpolant(R1',Z1',PHI1',mask') ;
clear mask
disp('interpolating...')
for iphi=1:length(phi)
mask(:,:,iphi)=Mask(R,Z,ones(size(Z)).*phi(iphi));
end
clear Mask
%%
if length(RAD)>1
disp('### Radiation')
disp('creating scatteredinterpolant...')
rad=scatteredInterpolant(Rp',Zp',PHIp',RADp') ;
clear RADp
disp('interpolating...')
for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldPrad(1,1,:,:,iphi)=rad(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
end
    clear rad
end

%%
disp('### Density')

for i=1:Ztot
    disp(['...density ',num2str(i),'/',num2str(size(Np,1))])
    disp('creating scatteredinterpolant...')
    %eval(['n',num2str(i),'=scatteredInterpolant(Rp'',Zp'',PHIp'',Np(',num2str(i),',:)'');']) ;
    n=scatteredInterpolant(Rp',Zp',PHIp',Np(i,:)') ;
    disp('interpolating...')
    for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldN(i,1,:,:,iphi)=n(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
    end
    clear n
end
clear Np
%%
disp('### Mach number')
disp('creating scatteredinterpolant...')
m=scatteredInterpolant(Rp',Zp',PHIp',Mp') ;
clear Mp
disp('interpolating...')
for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldM(1,1,:,:,iphi)=m(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
end
clear m

%%
disp('### T_e')
disp('creating scatteredinterpolant...')
te=scatteredInterpolant(Rp',Zp',PHIp',Tep') ;
clear Tep
disp('interpolating...')
for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldTe(1,1,:,:,iphi)=te(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
end
clear te
%%
disp('### T_i')
disp('creating scatteredinterpolant...')
ti=scatteredInterpolant(Rp',Zp',PHIp',Tip') ;
clear Tip
disp('interpolating...')
for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldTi(1,1,:,:,iphi)=ti(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
    end
    clear ti
%%
disp('### n_N')
disp('creating scatteredinterpolant...')
nn=scatteredInterpolant(Rpn',Zpn',PHIpn',Npn') ;
clear Npn
disp('interpolating...')
for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldNn(1,1,:,:,iphi)=nn(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
    end
    clear nn
%%
disp('### T_n')
disp('creating scatteredinterpolant...')
tn=scatteredInterpolant(Rpn',Zpn',PHIpn',Tpn') ;
clear Tpn
disp('interpolating...')
for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldTn(1,1,:,:,iphi)=tn(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
    end
    clear tn
%%
disp('### v_imp')

if Ztot>1
    for i=1:Ztot-1
        disp(['...impurity velocity ',num2str(i),'/',num2str(size(MIp,1))])
        disp('creating scatteredinterpolant...')
        %eval(['mi',num2str(i),'=scatteredInterpolant(Rp'',Zp'',PHIp'',MIp(',num2str(i),',:)'');']) ;
        mi=scatteredInterpolant(Rp',Zp',PHIp',MIp(i,:)') ;
        disp('interpolating...')
        for iphi=1:length(phi)
    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldMI(i,1,:,:,iphi)=mi(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
    end
    clear mi
    end
end
clear MIp

%%
T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time to this point ',char(s)])
%% radiation part
%     R1=repmat(R',1,1,length(phi));    Z1=repmat(Z',1,1,length(phi)); dphi1=ones(size(Z1));
%     phi1=ones(size(Z1));
%     for i=1:length(phi)
%         phi1(:,:,i)=phi(i);
%     end
%     dA=mean(diff(R(1,:))).*mean(diff(Z(:,1))).*mean(diff(phi)).*pi/180;
%     rad_p=ones(size(R,2),size(R,1),length(phi)).*NaN;
%     fieldPrad=squeeze(fieldPrad);
%  R1=R1(~isnan(fieldPrad)); Z1=Z1(~isnan(fieldPrad)); phi1=phi1(~isnan(fieldPrad));
%  fieldPrad1=fieldPrad(~isnan(fieldPrad)); 
% for k=1:length(phi)
%   k
%     dphi=phi1-phi(k);
%     cond1=dphi>round(max(phi))/2;
%     dphi(cond1>0)=dphi(cond1>0)-round(max(phi));
%     cond2=dphi<-round(max(phi))/2;
%     dphi(cond2>0)=dphi(cond2>0)+round(max(phi));
%     
%     
% for i=1:size(R,2)
%     for j=1:size(R,1)
%          D=(R1-R(j,i)).^2+(Z1-Z(j,i)).^2+(R1.*dphi.*pi/180).^2;%[m^2]
%          D(D==0)=NaN;
% rad_p(i,j,k)=nansum(fieldPrad1./(4.*pi.*D).*dA.*R1);%[W/m^2]
%     end 
% end
% end
% 


%%

mesh(1).xcells=R';
mesh(1).ycells=Z';
mesh(1).phicells=phi;
disp('saving....')
cd ..
save(fnameout,'mesh','fieldN','fieldTe','fieldTi','fieldNn','fieldTn','-v7.3')
disp(['saved in ',fnameout])

if Ztot>1
cs=9.79e3*sqrt(fieldTe);
for iz=1:Ztot-1
    fieldMI(iz,:,:,:,:)=fieldMI(iz,:,:,:,:)./cs(1,:,:,:,:);
end
save(fnameout,'fieldMI','-append')
disp('Impurity mach number added to file')
end

%% read kmag output
clear a
a=readmatrix(kmagfile,'FileType','text','NumHeaderLines',29);
a=a(:,1:3);
%%
% EMC3 always assumes that a positive Mach number represents the positive phi direction, or the counter-clockwise direction.
% However, KMAG gives the clockwise magnetic field
BR1=-reshape(a(:,1), size(R,1),size(R,2),length(phi));
Bphi1=-reshape(a(:,2), size(R,1),size(R,2),length(phi));
BZ1=-reshape(a(:,3), size(R,1),size(R,2),length(phi));
B1=sqrt(BR1.^2+BZ1.^2+Bphi1.^2);
clear B BR BZ Bphi
for i=1:size(B1,3)
    B(:,:,i)=squeeze(B1(:,:,i))';
    BR(:,:,i)=squeeze(BR1(:,:,i))';
    BZ(:,:,i)=squeeze(BZ1(:,:,i))';
    Bphi(:,:,i)=squeeze(Bphi1(:,:,i))';
end
clear BR1 BZ1 B1 Bphi1
%%
fieldmR(1,1,:,:,:)=squeeze(fieldM(1,1,:,:,:)).*BR./B;
fieldmZ(1,1,:,:,:)=squeeze(fieldM(1,1,:,:,:)).*BZ./B;
fieldmphi(1,1,:,:,:)=squeeze(fieldM(1,1,:,:,:)).*Bphi./B;
%%
save(fnameout,'fieldmR','fieldmZ','fieldmphi','-append','-v7.3')
disp(' M and B fields updated')
if Ztot>1
for i=1:Ztot-1
fieldmIR(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*BR./B;
fieldmIZ(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*BZ./B;
fieldmIphi(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*Bphi./B;
end
save(fnameout,'fieldmIR','fieldmIZ','fieldmIphi','-append')
disp('M field updated for impurities')
end

%% determine parameters for dust simulation
clear R1 R2
R1=2.743; R2=4.548;
kz=dsearchn(mesh(1).ycells(1,:)',0);
kr1=dsearchn(mesh(1).xcells(:,kz),R1);
kr2=dsearchn(mesh(1).xcells(:,kz),R2);
kphi=dsearchn(mesh(1).phicells',36/2);

n0=mean(fieldN(1,1,[kr1 kr2],kz,kphi));
T0=mean(fieldTe(1,1,[kr1 kr2],kz,kphi));
B0=2.75;
param.R0=3.6;
save(fnameout,'n0','T0','B0','param','-append','-v7.3')
T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time ',char(s)])
exit
catch e
    try disp(e.stack); catch; end
    disp(e.identifier)
    disp(e.message)
    disp('exiting')
    T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time to this point ',char(s)])
  %   exit; 
end

%% Explanation from Yuhe Feng
% The geometric cells
% range from 0 to MESH_P_OS(NZONET)-1, where NZONET is the total number of
% mesh blocks/zones and MESH_P_OS(NZONET) is the total geometric cell
% number.
% The 1D address of a geometric cell, IG, is linked to its address
% in real space, which is labelled by IZ (mesh block nr.), I (radial
% index), J (poloidal index) and k (toroidal), uniquely by
% IG=I+J*ZON_RADI(IZ)+K*ZON_RADI(IZ)*ZON_POLO(IZ)+MESH_P_OS(IZ)

% where ZON_RADI(IZ) and ZON_POLO(IZ) are respectively the total radial
% and poloidal mesh points of a standard mesh block IZ and MESH_P_OS(IZ)
% stands for the offset of the first cell number in this mesh block.
% A physical cell IC is linked to a gemometric cell IG by

% IC = IDCELL(IG)
% where the 1D array IDCELL is stored in the file CELL_GEO.
% Each gemotric cell [IZ,I,J,K] is formed by eight grid points. The
% addresses of the corner points are in the 1D arrays of RG(0:6421339) and
% ZG(0:6421339) are respectively
% I1= I  + J  *SRF_RADI(IZ)+K*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I2= I+1+ J  *SRF_RADI(IZ)+K*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I3= I+1+ J+1*SRF_RADI(IZ)+K*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I4= I  + J+1*SRF_RADI(IZ)+K*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I5= I  + J  *SRF_RADI(IZ)+(K+1)*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I6= I+1+ J  *SRF_RADI(IZ)+(K+1)*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I7= I+1+ J+1*SRF_RADI(IZ)+(K+1)*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% I8= I  + J+1*SRF_RADI(IZ)+(K+1)*SRF_RADI(IZ)*SRF_POLO(IZ)+GRID_P_OS(IZ)
% where SRF_RADI(IZ)=ZON_RADI(IZ)+1, SRF_POLO(IZ)=ZON_POLO(IZ)+1 and
% GRID_P_OS(IZ) is the offset of grid points, in analogous to MESH_P_OS
% for geometric cells. The eight points, together its phi angle, define
% its location in real space.

