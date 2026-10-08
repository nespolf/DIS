clear all

try
tic;
test=1;
fnameout=['DIS_pbg_EMC3.mat'];

if test

cd ~/EMC3_post/SRO_Lmode_20deg_shapedWall/run/EMC3_OUTPUT/
inputgeofile='/Users/fnespoli/EMC3_post/SRO_Lmode_20deg_shapedWall/geometry/input.geo';
gridfile='/Users/fnespoli/EMC3_post/SRO_Lmode_20deg_shapedWall/geometry/grid3D.dat';

cd ~/EMC3_post/Aug27_Ne1.2_Florian/Ne1.2_0826_6/EMC3_OUTPUT/
inputgeofile='/Users/fnespoli/EMC3_post/Aug27_Ne1.2_Florian/geometry/input.geo';
gridfile='/Users/fnespoli/EMC3_post/Aug27_Ne1.2_Florian/geometry/grid3D.dat';

cd ~/EMC3_post/Limiter_Case/run_new/EMC3_OUTPUT/
inputgeofile='/Users/fnespoli/EMC3_post/Limiter_Case/geometry/input.geo';
gridfile='/Users/fnespoli/EMC3_post/Limiter_Case/geometry/grid3D.dat';


else
    curdir=pwd;
inputgeofile=[curdir,'/input/input.geo'];
gridfile=[curdir,'/','GEO/G']; 
fnameout=[curdir,'/',fnameout];
end



%% read grid info from input.geo
disp('read grid info...')
fid=fopen(inputgeofile);
for i=1:7
    a=fgetl(fid);
end

NZ=str2num(a); IZ=1:NZ;
clear Nr Nth Nphi
a=fgetl(fid);
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

clear R Z PHI Phi
fid=fopen(gridfile);
 R1=[]; Z1=[];  step=1;  Nzones=11;

    for iz=1:NZ
    
    
        n=str2num(fgetl(fid));
        for iphi=1:n(3)
            R1=[]; Z1=[];
        a=str2num(fgetl(fid));
        Phi(iphi,iz)=a; 

        for i=1:ceil(n(1)*n(2)/4)
        a=str2num(fgetl(fid));
        R1=[R1 a]; 
        end
        for i=1:ceil(n(1)*n(2)/4)
        a=str2num(fgetl(fid));
        Z1=[Z1 a]; 
        end
        R(iz).d(iphi,:)=R1;
        Z(iz).d(iphi,:)=Z1;
        PHI(iz).d(iphi,:)=ones(size(R1)).*Phi(iphi,iz);
        end
        
    end

   
fclose(fid);

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

%%
figure

for i=1%:size(R(iz).d,1)
    for iz=1:length(R)
    scatter(R(iz).d(i,:),Z(iz).d(i,:)); axis equal; hold on
    end
    hold off
     pause(0.1)
end
%% magnetic field components from the grid
figure
RG20=[]; ZG20=[]; PHIG20=[]; BR20=[]; BZ20=[]; BPHI20=[];
for iz=1:length(Nr)
     a=R(iz).d';
    b=Z(iz).d';
    c=PHI(iz).d';
RG1=reshape(a,Nr(iz),Nth(iz),Nphi(iz));
ZG1=reshape(b,Nr(iz),Nth(iz),Nphi(iz));
PHIG1=reshape(deg2rad(c),Nr(iz),Nth(iz),Nphi(iz));
dphi=unique(diff(PHIG)); dphi=dphi(2);
clear dR dZ dL dPHI RG2 ZG2 PHIG2
for k=1:Nphi(iz)-1
    dR(:,:,k)=RG1(:,:,k+1)-RG1(:,:,k);
    dZ(:,:,k)=ZG1(:,:,k+1)-ZG1(:,:,k);
    dPHI(:,:,k)=0.5.*(RG1(:,:,k+1)+RG1(:,:,k)).*(PHIG1(:,:,k+1)-PHIG1(:,:,k));

    RG2(:,:,k)=0.5.*(RG1(:,:,k+1)+RG1(:,:,k));
    ZG2(:,:,k)=0.5.*(ZG1(:,:,k+1)+ZG1(:,:,k));
    PHIG2(:,:,k)=0.5.*(PHIG1(:,:,k+1)+PHIG1(:,:,k));
end
 dL=sqrt(dR.^2+dZ.^2+dPHI.^2);
 BR=dR./dL; 
 BZ=dZ./dL;
 BPHI=dPHI./dL;
 RG20=[RG20 RG2(:)'];
 ZG20=[ZG20 ZG2(:)'];
 PHIG20=[PHIG20 PHIG2(:)'];
 BR20=[BR20 BR(:)'];
 BZ20=[BZ20 BZ(:)'];
 BPHI20=[BPHI20 BPHI(:)'];

%figure; 
for i=1%:size(BPHI,3)
    subplot(1,3,1)
    pcolor(squeeze(RG2(:,:,i)),squeeze(ZG2(:,:,i)),squeeze(BR(:,:,i))); shading flat; axis equal;hold on
    title('B_R'); set(gca,'fontsize',14); colorbar
    subplot(1,3,2)
    pcolor(squeeze(RG2(:,:,i)),squeeze(ZG2(:,:,i)),squeeze(BZ(:,:,i))); shading flat; axis equal;hold on
    title('B_Z'); set(gca,'fontsize',14); colorbar
    subplot(1,3,3)
    pcolor(squeeze(RG2(:,:,i)),squeeze(ZG2(:,:,i)),squeeze(BPHI(:,:,i))); shading flat; axis equal; hold on
    title('B_{PHI}'); set(gca,'fontsize',14); colorbar
    colormap jet
    pause(0.1)
end


end
RG2=RG20; ZG2=ZG20; PHIG2=PHIG20; BR=BR20; BZ=BZ20; BPHI=BPHI20;
clear RG20 ZG20 PHIG20 BPHI20 BR20 BZ20
%% read cell geo
disp('read cell geo...')

clear IDCELL
out=readmatrix('../../geometry/IDCELL','FileType','text','NumHeaderLines',1)';
IDCELL=out(:);
IDCELL=IDCELL(~isnan(IDCELL));
fid=fopen('../../geometry/IDCELL');
Ncell=str2num(fgetl(fid));
if isempty(IDCELL)%isnan(IDCELL) case from stellar with CELL_GEO on a single line
disp('a')
IDCELL=str2num(fgetl(fid));
end
fclose(fid);
%% read lg-cell
% disp('read cell length...')
% 
% clear LGCELL
% out=readmatrix('../LG_CELL','FileType','text','NumHeaderLines',0)';
% 
% out1=out;
% LGCELL=out1(:);
% LGCELL=LGCELL(~isnan(LGCELL));
% 
% 
% if ~test cd OUTPUT; end
%% read Temperature
%  %cd /Users/fnespoli/EMC3_post/W7-X_Kawamura/data-360
% 
% 
% disp('read Temperature...')
% 
% clear out
% out=readmatrix('TE_TI','FileType','text')';
% out=out(:);
% idx=find(isnan(out));
% out=out(isnan(out)<1);
% if isempty(idx) idx=length(out)/2+1; end
% TE=out(1:idx(1)-1);
% TI=out(idx(1):end);
%% read Temperature 2nd method
 %cd /Users/fnespoli/EMC3_post/W7-X_Kawamura/data-360
disp('read Temperature...')
clear out
out=readmatrix('TE_TI','FileType','text')';
out=out(:);
out=out(isnan(out)<1);
if length(out)/2==Ncell(2)
    TE=out(1:Ncell(2));
    TI=out(Ncell(2)+1:end);
else
    disp('this method is wrong')
end
%% read density
% disp('read density...')
% 
% clear out N
% out=readmatrix('DENSITY','FileType','text','NumHeaderLines',0)';
% 
% out=out(:);
% idx=find(isnan(out));
% 
% idxprev=0; k=0; k1=0; idx1(1)=0;
% for i=1:length(idx)
%     if idx(i)-idxprev>Ncell(2)
%         k=k+1;
%         idx1(k)=idx(i);
%         idxprev=idx(i);
%     end
%     if i>1 && idx(i)-idx(i-1)>Ncell(2)
%         k1=k1+1;
%         idx2(k)=idx(i-1);
%     end
% end
% if ~exist('idx2') idx2=0; end
% for i=1:length(idx1)
%     N(i,:)=out(idx2(i)+1:idx1(i)-1);
% end
% 
% if isempty(N) N=out'; end
%% read density 2nd method
try
disp('read density...')
clear out N
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
if size(N,2)<Ncell(2) a=N(i,Ncell(2)); end
catch
disp('did not work')
     disp('read Density 3rd method...')
clear out N 
out=readmatrix('DENSITY','FileType','text','NumHeaderLines',0)';
Ztot=floor(length(out(:))/Ncell(2));
Nlines=Ncell(2)/6;%for each species
if ~isinteger(Nlines) Nlines=ceil(Nlines); end
for i=1:Ztot
A=out(:,(i-1)*Nlines+1+i:i*Nlines+i);
A=A(:);
if Nlines*6>Ncell(2) A=A(1:Ncell(2)); end
N(i,:)=A;
end
end
%% read mach number
disp('read mach number...')

clear out M
out=readmatrix('MACH_NUMBER','FileType','text')';
out=out(:);
out(isnan(out)==1)=0; %fix for when mach number is too small e.g. E-270 formatting gets wrong and E disappears, gets read as NaN 
M=out(1:Ncell(2));

%% read neutral density
disp('read neutral density...')
try
clear out
out=readmatrix('../EIRENE_OUTPUT/DENSITY_A','FileType','text')';
out=out(:);
out=out(isnan(out)<1);
Nn=out(:);
catch
    Nn=[]; disp('no neutral density')
end
%% read neutral Temperature
disp('read neutral Temperature...')
try
clear out
out=readmatrix('../EIRENE_OUTPUT/TEMPERATURE_A','FileType','text')';
out=out(:);
out=out(isnan(out)<1);
Tn=out(:);
catch
    Tn=[]; disp('no neutral temperature')
end
%% read impurity velocity
disp('read impurity velocity...')

%if size(N,1)>1
 try  
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
 catch
     disp('no impurity velocity')
    MI=0;
end



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
    if IC<=Ncell(2) && IC>0
        Np(:,k)=N(:,IC);
       % MIp(:,k)=MI(:,IC);
        Mp(k)=M(IC);
        Tep(k)=TE(IC);
        Tip(k)=TI(IC);
 %       RADp(k)=RAD(IC);
        Rp(k)=R1(i);  Zp(k)=Z1(i); PHIp(k)=PHI1(i); mask(i)=1;
        k=k+1;
    else
        R1(i)=mean(RG(idx));
        Z1(i)=mean(ZG(idx)); PHI1(i)=mean(PHIG(idx));  mask(i)=NaN;
    end
    if IC<=Ncell(3) && ~isempty(Nn)
        Npn(k1)=Nn(IC); 
        if ~isempty(Tn) Tpn(k1)=Tn(IC);end
        Rpn(k1)=R1(i);  Zpn(k1)=Z1(i); PHIpn(k1)=PHI1(i);
        k1=k1+1;
    else
    end
end


%%
% 
cond=PHIp==(PHIp(dsearchn(PHIp',0)));
figure;
scatter(Rp(cond>0),Zp(cond>0),5,(Np(1,cond>0)),'fill'); axis equal

%%

clear R Z RG ZG PHIG N TE TI M MI I1 I2 I3 I4 I5 I6 I7 I8  I11 I21 I31 I41 I51 I61 I71 I81 IG IDCELL Nn Tn
Rp=Rp/100; Zp=Zp/100; R1=R1/100; Z1=Z1/100; Np=Np.*1e6;

Ztot=size(Np,1);

Nr=250; Nth=500;
%Nr=250; Nth=250;
%Nr=99; Nth=101;

R=linspace(min(Rp),max(Rp),Nr);
Z=linspace(-max(abs(Zp)),max(abs(Zp)),Nth);
phi=linspace(min(PHIp),max(PHIp),32);
phiu=unique(PHI1); phi=phiu;

phi=linspace(-20,20,32);% for FT grid
phiu=unique(PHI1); 
cond=(phiu>-25).*(phiu<25);
phiu=phiu(cond>0);

[R Z]=meshgrid(R,Z);

clear fieldN fieldTe fieldTi fieldM fieldNn fieldTn fieldMI fieldPrad
if exist('Npn')
Rpn=Rpn/100; Zpn=Zpn/100; Npn=Npn.*1e6;

A=Rpn(PHIpn==min(PHIpn)); B=Zpn(PHIpn==min(PHIpn));
k=boundary(A(:),B(:));
A=A(k); B=B(k);
disp('interpolating mask')

maskn=double(inpolygon(R,Z,A,B));
maskn(maskn==0)=NaN;
end
%%
disp('### Mask')
disp('creating scatteredinterpolant, interpolating...')
for k=1:length(phiu)
    cond=PHI1==phiu(k);
    Mask=scatteredInterpolant(R1(cond>0)',Z1(cond>0)',mask(cond>0)') ;
    mask1(:,:,k)=Mask(R,Z);
end
mask=mask1;
clear Mask
% %%
% if length(RAD)>1
% disp('### Radiation')
% disp('creating scatteredinterpolant...')
% rad=scatteredInterpolant(Rp',Zp',PHIp',RADp') ;
% clear RADp
% disp('interpolating...')
% for iphi=1:length(phi)
%     %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
%     fieldPrad(1,1,:,:,iphi)=rad(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
% end
%     clear rad
% end
% 
%%
disp('### Density')

for i=1%:Ztot
    disp(['...density ',num2str(i),'/',num2str(size(Np,1))])
    disp('creating scatteredinterpolant, interpolating...')
    for k=1:length(phiu)
    cond=PHIp==phiu(k);
    n=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Np(i,cond>0)') ;
    fieldN(i,1,:,:,k)=n(R,Z)'.*mask(:,:,k)';
    clear n
    end
end
clear Np
%%
disp('### Mach number')

disp('creating scatteredinterpolant, interpolating...')
    for k=1:length(phiu)
    cond=PHIp==phiu(k);
    m=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Mp(1,cond>0)') ;
    fieldM(1,1,:,:,k)=m(R,Z)'.*mask(:,:,k)';
    clear m
    end

clear Mp

%%
disp('### T_e')
disp('creating scatteredinterpolant, interpolating...')
    for k=1:length(phiu)
    cond=PHIp==phiu(k);
    te=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Tep(1,cond>0)') ;
    fieldTe(1,1,:,:,k)=te(R,Z)'.*mask(:,:,k)';
    clear te
    end
clear te
%%
disp('### T_i')
disp('creating scatteredinterpolant, interpolating...')
    for k=1:length(phiu)
    cond=PHIp==phiu(k);
    ti=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Tip(1,cond>0)') ;
    fieldTi(1,1,:,:,k)=ti(R,Z)'.*mask(:,:,k)';
    clear ti
    end
clear Tip

%%
if exist('Npn')
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
end
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
% clear MIp

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

cs=9.79e3*sqrt(fieldTe);
if ~exist('fieldMI') fieldMI=[]; else
for iz=1:Ztot-1
    fieldMI(iz,:,:,:,:)=fieldMI(iz,:,:,:,:)./cs(1,:,:,:,:);
end
end
mesh(1).xcells=R';
mesh(1).ycells=Z';
%mesh(1).phicells=phi;
mesh(1).phicells=phiu;
disp('saving....')
if ~exist('fieldMI') fieldMI=[]; end
if ~exist('fieldNn') fieldNn=[]; fieldTn=[];end
save(fnameout,'mesh','fieldN','fieldTe','fieldTi','fieldM','fieldMI','fieldNn','fieldTn','-v7.3')
disp(['saved in ',fnameout])


%% interpolate B field
warning off
disp('### B')
disp('creating scatteredinterpolant...')
phiu1=unique(PHIG2(:)); 
cond=(phiu1>deg2rad(-25)).*(phiu1<deg2rad(25));%downsizing FT grid
phiu1=phiu1(cond>0);
if ~phiu==phiu1 disp('problem here'); end
RG22=RG2(:); ZG22=ZG2(:);  PHIG22=PHIG2(:); BR2=BR(:); BZ2=BZ(:); BPHI2=BPHI(:);

for k=1:length(phiu1)
    cond=PHIG22==phiu1(k);
br=scatteredInterpolant(RG22(cond>0)./100,ZG22(cond>0)./100,BR2(cond>0)) ;
bz=scatteredInterpolant(RG22(cond>0)./100,ZG22(cond>0)./100,BZ2(cond>0)) ;
bphi=scatteredInterpolant(RG22(cond>0)./100,ZG22(cond>0)./100,BPHI2(cond>0)) ;

    %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
    fieldBR(1,1,:,:,k)=br(R,Z)'.*mask(:,:,k)';
    fieldBZ(1,1,:,:,k)=bz(R,Z)'.*mask(:,:,k)';
    fieldBphi(1,1,:,:,k)=bphi(R,Z)'.*mask(:,:,k)';
end
    
disp('hop')
%% from half period to full period
if min(phi)>0 %half period
    disp('from half period to full period')
load(fnameout)
phi1=[-phi(end:-1:1) phi];
for k=1:length(phi)
 fieldN1(:,:,:,:,k)=fieldN(:,:,:,end:-1:1,end-(k-1));
 fieldN1(:,:,:,:,length(phi)+k)=fieldN(:,:,:,:,k);

 fieldTe1(:,:,:,:,k)=fieldTe(:,:,:,end:-1:1,end-(k-1));
 fieldTe1(:,:,:,:,length(phi)+k)=fieldTe(:,:,:,:,k);

 fieldTi1(:,:,:,:,k)=fieldTi(:,:,:,end:-1:1,end-(k-1));
 fieldTi1(:,:,:,:,length(phi)+k)=fieldTi(:,:,:,:,k);

 fieldM1(:,:,:,:,k)=-fieldM(:,:,:,end:-1:1,end-(k-1));
 fieldM1(:,:,:,:,length(phi)+k)=fieldM(:,:,:,:,k);

 fieldBR1(:,:,:,:,k)=fieldBR(:,:,:,end:-1:1,end-(k-1));
 fieldBR1(:,:,:,:,length(phi)+k)=fieldBR(:,:,:,:,k);

 fieldBZ1(:,:,:,:,k)=fieldBZ(:,:,:,end:-1:1,end-(k-1));
 fieldBZ1(:,:,:,:,length(phi)+k)=fieldBZ(:,:,:,:,k);

 fieldBphi1(:,:,:,:,k)=fieldBphi(:,:,:,end:-1:1,end-(k-1));
 fieldBphi1(:,:,:,:,length(phi)+k)=fieldBphi(:,:,:,:,k);
end

fieldN=fieldN1;
fieldTe=fieldTe1;
fieldTi=fieldTi1;
fieldM=fieldM1;
fieldBR=fieldBR1;
fieldBZ=fieldBZ1;
fieldBphi=fieldBphi1;

clear fieldBphi1 fieldBR1 fieldBZ1 fieldN1 fieldM1 fieldTe1 fieldTi1
%%
if test
figure;
for k=1:size(fieldN,5)
   % pcolor(R,Z,squeeze(mask(:,:,k))); shading flat
    pcolor(R,Z,squeeze(fieldN(1,1,:,:,k))'); shading flat

    colormap jet; caxis([-0.5 0.5])
    pause
end
end
%%

mesh.phicells=phi1; 
save(fnameout,'mesh','fieldN','fieldTe','fieldTi','fieldM','fieldMI','fieldNn','fieldTn','-append')
disp(['fields updated in ',fnameout])
end
disp('hop')
%%
clear BR BZ PBHI br bz bphi
fieldB=sqrt(fieldBR.^2+fieldBZ.^2+fieldBphi.^2);

%%
fieldmR=fieldM.*fieldBR./fieldB;
fieldmZ=fieldM.*fieldBZ./fieldB;
fieldmphi=fieldM.*fieldBphi./fieldB;
fieldmIR=[]; fieldmIZ=[]; fieldmIphi=[];
% 
% for i=1:Ztot-1
% fieldmIR(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*BR./B;
% fieldmIZ(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*BZ./B;
% fieldmIphi(i,1,:,:,:)=squeeze(fieldMI(i,1,:,:,:)).*Bphi./B;
% end




%%
save(fnameout,'fieldmR','fieldmZ','fieldmphi','fieldmIR','fieldmIZ','fieldmIphi','-append','-v7.3')
disp(' M and B fields updated')


%% determine parameters for dust simulation
clear R1 R2
% R1=2.743; R2=4.548;
% kz=dsearchn(mesh(1).ycells(1,:)',0);
% kr1=dsearchn(mesh(1).xcells(:,kz),R1);
% kr2=dsearchn(mesh(1).xcells(:,kz),R2);
% kphi=dsearchn(mesh(1).phicells',36/2);
% 
% n0=mean(fieldN(1,1,[kr1 kr2],kz,kphi));
% T0=mean(fieldTe(1,1,[kr1 kr2],kz,kphi));
n0=nanmean(nanmean(nanmean(squeeze(fieldN(1,:,:,:,:)))));
T0=nanmean(nanmean(nanmean(squeeze(fieldTe))));
B0=3.5;
param.R0=1.8;
save(fnameout,'n0','T0','B0','param','-append','-v7.3')
T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time ',char(s)])
if ~test  exit; end 
catch e
    try disp(e.stack); catch; end
    disp(e.identifier)
    disp(e.message)
    disp(['file', e.stack.file, ' line ', num2str(e.stack.line)])
    disp('exiting')
    T1=toc;
s = seconds(T1);
s.Format = 'hh:mm:ss';
disp(['total time to this point ',char(s)])
  if ~test  exit; end 
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

