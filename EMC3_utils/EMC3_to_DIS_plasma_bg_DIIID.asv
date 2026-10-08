clear all
try
    curdir=pwd;
tic;
fnameout=['/Users/fnespoli/scripts/DIS/plasma_background/DIS_pbg_EMC3_DIIID_USN.mat'];
% inputgeofile='/Users/fnespoli/EMC3_post/DIIID/input.geo';
% gridfile='/Users/fnespoli/EMC3_post/DIIID/geometry/grid3D.dat';
% idcellfile='/Users/fnespoli/EMC3_post/DIIID/geometry/IDCELL';
% emc3resultsfolder='/Users/fnespoli/EMC3_post/DIIID/output/';
% kmagfile='/p/emc3-eirene/fnespoli/kmag/output.txt';

inputgeofile='/Users/fnespoli/EMC3_post/DIIID/179894_3500ms/EMC3_GRID/input.geo';
gridfile='/Users/fnespoli/EMC3_post/DIIID/179894_3500ms/EMC3_GRID/fullgridwall.3d';
idcellfile='/Users/fnespoli/EMC3_post/DIIID/179894_3500ms/EMC3_GRID/CELL_GEO';
emc3resultsfolder='/Users/fnespoli/EMC3_post/DIIID/179894_3500ms/OUTPUT/';
kmagfile='/Users/fnespoli/EMC3_post/DIIID/179894_3500ms/Bfield_sampled/Bfield.dat';
kmaggridfile='/Users/fnespoli/EMC3_post/DIIID/179894_3500ms/Bfield_sampled/grid.dat';

inputgeofile='/Users/fnespoli/EMC3_post/DIIID/176750_2950ms_APS/NE15E18P6MW/run/input.geo';
gridfile='/Users/fnespoli/EMC3_post/DIIID/176750_2950ms_APS/geometry/grid3D.dat';
idcellfile='/Users/fnespoli/EMC3_post/DIIID/176750_2950ms_APS/geometry/CELL_GEO';
emc3resultsfolder='/Users/fnespoli/EMC3_post/DIIID/176750_2950ms_APS/NE15E18P6MW/run/OUTPUT/';
kmagfile='/Users/fnespoli/EMC3_post/DIIID/176750_2950ms_APS/Bfield_five_columns.dat';
kmaggridfile='/Users/fnespoli/EMC3_post/DIIID/176750_2950ms_APS/sample_grid_used_for_Bfield_sampling.dat';

%% read grid info from input.geo
disp('read grid info...')
fid=fopen(inputgeofile);
for i=1:3
    a=fgetl(fid);
    disp(a)
end
NZ=str2num(a); IZ=1:NZ;
clear Nr Nth Nphi
for i=1:NZ
    a=fgetl(fid);     
     k=findstr(a,'!'); %for LSN
     if ~isempty(k)
     a=a(1:k-1);
     end
    a=str2num(a);
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

clear R Z PHI Phi
fid=fopen(gridfile,'r');
A=fscanf(fid,'%f');
fclose(fid);
j=0; iz=1; i=1;
for iz=1:NZ
    j=j+3;
    for i=1:Nphi(iz)
    j=j+1;
    Phi(i,iz)=A(j);
    R1=A(j+1:j+Nr(iz).*Nth(iz));
    Z1=A(j+Nr(iz).*Nth(iz)+1:j+2.*Nr(iz).*Nth(iz));
    R(iz).d(i,:)=R1;
    Z(iz).d(i,:)=Z1;
    PHI(iz).d(i,:)=ones(size(R1)).*Phi(i,iz);
    j=j+2.*Nr(iz).*Nth(iz);
    end
end

figure;
for iz=1:NZ
    scatter(R(iz).d(1,:),Z(iz).d(1,:)); axis equal; hold on
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
% out=readmatrix(idcellfile,'FileType','text','NumHeaderLines',1)';
% out1=out;
% IDCELL=out1(:)+1;
fid=fopen(idcellfile);
Ncell=str2num(fgetl(fid));
IDCELL=str2num(fgetl(fid));
fclose(fid);

%% read Temperature
cd(emc3resultsfolder)

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

clear out N
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
if ~exist('idx2') idx2=0; end
for i=1:length(idx1)
    N(i,:)=out(idx2(i)+1:idx1(i)-1);
end
%% read neutral density
disp('read neutral density...')

clear out
out=readmatrix('DENSITY_A','FileType','text')';
out=out(:);
out=out(isnan(out)<1);
Nn=out(:);
% %% read neutral Temperature
% disp('read neutral Temperature...')
% 
% clear out
% out=readmatrix('TEMPERATURE_A','FileType','text')';
% out=out(:);
% out=out(isnan(out)<1);
% Tn=out(:);
% 

%% read mach number
disp('read mach number...')

clear out M
out=readmatrix('MACH_NUMBER','FileType','text')';
out=out(isnan(out)<1);
M=out(:);

%% read impurity velocity
% disp('read impurity velocity...')
% 
% if size(N,1)>1
%     clear out MI idx idx1 idx2
%     out=readmatrix('IMPURITY_VELOCITY_2','FileType','text')';
%     
%     out=out(:);
%     idx=find(isnan(out));
%     idxprev=0; k=0; k1=0; idx1(1)=0;
%     for i=1:length(idx)
%         if idx(i)-idxprev>Ncell(2)
%             k=k+1;
%             idx1(k)=idx(i);
%             idxprev=idx(i);
%         end
%         if i>1 && idx(i)-idx(i-1)>Ncell(2)
%             k1=k1+1;
%             idx2(k)=idx(i-1);
%         end
%     end
%     
%     for i=1:length(idx1)
%         MI(i,:)=out(idx2(i)+1:idx1(i)-1);
%     end
% end

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
k=1; k1=1;
for i=1:length(IG)%Ncell(2)%length(IC)
    IC=int32(IDCELL(IG(i)));
    idx=[I11(i) I21(i) I31(i) I41(i) I51(i) I61(i) I71(i) I81(i)];
    R1(i)=mean(RG(idx));
    Z1(i)=mean(ZG(idx));
    PHI1(i)=mean(PHIG(idx));
    if IC<=Ncell(2)
        Np(:,k)=N(:,IC);
    %    MIp(:,k)=MI(:,IC);
        Mp(k)=M(IC);
        Tep(k)=TE(IC);
        Tip(k)=TI(IC);
        if length(RAD)==1 RADp(k)=0;
        else RADp(k)=RAD(IC);
        end
        Rp(k)=R1(i);  Zp(k)=Z1(i); PHIp(k)=PHI1(i); mask(i)=1;
        k=k+1;
    else
        R1(i)=mean(RG(idx));
        Z1(i)=mean(ZG(idx)); PHI1(i)=mean(PHIG(idx));  mask(i)=NaN;
    end
    if IC<=Ncell(3)
        Npn(k1)=Nn(IC);% Tpn(k1)=Tn(IC);
        Rpn(k1)=R1(i);  Zpn(k1)=Z1(i); PHIpn(k1)=PHI1(i);
        k1=k1+1;
    else
    end
end


%%
clear R Z RG ZG PHIG N TE TI M MI I1 I2 I3 I4 I5 I6 I7 I8  I11 I21 I31 I41 I51 I61 I71 I81 IG IDCELL Nn Tn
Rp=Rp/100; Zp=Zp/100; R1=R1/100; Z1=Z1/100; Np=Np.*1e6;
Rpn=Rpn/100; Zpn=Zpn/100; Npn=Npn.*1e6;
Ztot=size(Np,1);

Nr=200; Nth=370;
Nr=200; Nth=400;


R=linspace(min(Rp),max(Rp),Nr);
Z=linspace(min(Zp),max(Zp),Nth);
phi=linspace(min(PHIp),max(PHIp),64);
phi=0;

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
% disp('### T_n')
% disp('creating scatteredinterpolant...')
% tn=scatteredInterpolant(Rpn',Zpn',PHIpn',Tpn') ;
% clear Tpn
% disp('interpolating...')
% for iphi=1:length(phi)
%     %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
%     fieldTn(1,1,:,:,iphi)=tn(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
%     end
%     clear tn
%%
% disp('### v_imp')
% 
% if Ztot>1
%     for i=1:Ztot-1
%         disp(['...impurity velocity ',num2str(i),'/',num2str(size(MIp,1))])
%         disp('creating scatteredinterpolant...')
%         %eval(['mi',num2str(i),'=scatteredInterpolant(Rp'',Zp'',PHIp'',MIp(',num2str(i),',:)'');']) ;
%         mi=scatteredInterpolant(Rp',Zp',PHIp',MIp(i,:)') ;
%         disp('interpolating...')
%         for iphi=1:length(phi)
%     %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
%     fieldMI(i,1,:,:,iphi)=mi(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
%     end
%     clear mi
%     end
% end
% clear MIp

%%
T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time to this point ',char(s)])
% %% radiation part
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



%%

cs=9.79e3*sqrt(fieldTe);
% for iz=1:Ztot-1
%     fieldMI(iz,:,:,:,:)=fieldMI(iz,:,:,:,:)./cs(1,:,:,:,:);
% end

mesh(1).xcells=R';
mesh(1).ycells=Z';
mesh(1).phicells=phi;
disp('saving....')
cd ..
save(fnameout,'mesh','fieldN','fieldM','fieldTe','fieldTi','fieldNn','-v7.3')
disp(['saved in ',fnameout])


%% read kmag output
clear a B1 BR1 BZ1 Bphi1
%a=readmatrix(kmagfile,'FileType','text','NumHeaderLines',8);%LSN file
%a=readmatrix(kmagfile,'FileType','text','NumHeaderLines',29);%USN file
a=readmatrix(kmagfile,'FileType','text','CommentStyle',{'#'});
B1=a(:,1); BR1=a(:,2); BZ1=a(:,3); Bphi1=a(:,4);
clear a
%a=readmatrix(kmaggridfile,'FileType','text','NumHeaderLines',3);
a=readmatrix(kmaggridfile,'FileType','text','CommentStyle',{'#'});
%R1=a(:,1); Z1=a(:,2);
Nr1=length(unique(a(:,1))); Nth1=length(unique(a(:,2))); 
%R1=reshape(R1,Nth1,Nr1)';
%Z1=reshape(Z1,Nth1,Nr1)';
[R1 Z1]=meshgrid(unique(a(:,1)), unique(a(:,2)));
B1=reshape(B1,Nth1,Nr1)';
BR1=reshape(BR1,Nth1,Nr1)';
BZ1=reshape(BZ1,Nth1,Nr1)';
Bphi1=reshape(Bphi1,Nth1,Nr1)';
B=interp2(R1./100,Z1./100,B1',R,Z);
BR=interp2(R1./100,Z1./100,BR1',R,Z);
BZ=interp2(R1./100,Z1./100,BZ1',R,Z);
Bphi=interp2(R1./100,Z1./100,Bphi1',R,Z);
B=repmat(B',1,1,length(mesh(1).phicells));
BR=repmat(BR',1,1,length(mesh(1).phicells));
BZ=repmat(BZ',1,1,length(mesh(1).phicells));
Bphi=repmat(Bphi',1,1,length(mesh(1).phicells));

%%
% EMC3 always assumes that a positive Mach number represents the positive phi direction, or the counter-clockwise direction.
%%

fieldmR(1,1,:,:,:)=squeeze(fieldM(1,1,:,:,:)).*BR./B;
fieldmZ(1,1,:,:,:)=squeeze(fieldM(1,1,:,:,:)).*BZ./B;
fieldmphi(1,1,:,:,:)=squeeze(fieldM(1,1,:,:,:)).*Bphi./B;
% for i=1:Ztot-1
% fieldmIR(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*BR./B;
% fieldmIZ(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*BZ./B;
% fieldmIphi(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*Bphi./B;
% end
%%
save(fnameout,'fieldmR','fieldmZ','fieldmphi','-append','-v7.3')
disp(' M and B fields updated')
%% determine parameters for dust simulation
clear R1 R2
R1=1.123; R2=2.246;
kz=dsearchn(mesh(1).ycells(1,:)',0);
kr1=dsearchn(mesh(1).xcells(:,kz),R1);
kr2=dsearchn(mesh(1).xcells(:,kz),R2);
kphi=dsearchn(mesh(1).phicells',0);

n0=mean(fieldN(1,1,[kr1 kr2],kz,kphi));
T0=mean(fieldTe(1,1,[kr1 kr2],kz,kphi));
B0=1.984;
param.R0=1.721;
B0=2;  param.R0=1.74; n0=3e19; T0=250; %USN case
save(fnameout,'n0','T0','B0','param','-append','-v7.3')
T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time ',char(s)])
cd(curdir)
%exit
catch e
    try disp(e.stack); catch; end
    disp(e.identifier)
    disp(e.message)
    disp('exiting')
    T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time to this point ',char(s)])
cd(curdir)
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

