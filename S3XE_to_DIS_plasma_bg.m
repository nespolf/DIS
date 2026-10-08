clear all

try
tic;
test=1;
fnameout=['DIS_pbg_S3XE.mat'];

if test
addpath ~/scripts/
folder_name={'/Users/fnespoli/S3XE/WEST_D'};
 foldtot=length(folder_name);
 fnameout=[folder_name{1},'/',fnameout];
else
    curdir=pwd;
inputgeofile=[curdir,'/input/input.geo'];
gridfile=[curdir,'/','GEO/G']; 
fnameout=[curdir,'/',fnameout];
end


for  i=1:foldtot
 folder_list1{i}=[folder_name{i},'/run_dir/'];
end
param = s3xe_load_param(folder_list1{1});
[param,mesh] = s3xe_load_mesh(folder_list1{1},param,1);

%%
Ntload_force = -1; % Set to -1 to let matlab load all the available time steps
    plotbool =0;
    dtmin = 0;
    Nt0=0;dtload=1;
    fields_list={'n','T','neutrals/Sn','ExB/phi','G'};
    [param,mesh,metric,Bfield,tt,fields] = s3xe_load_all(folder_list1,fields_list,Nt0,dtload,plotbool,Ntload_force);

for iz=1:length(fields)
 fields(iz).M=fields(iz).G./fields(iz).n./sqrt((fields(iz).T(1,:,:,:,:)+fields(iz).T(2,:,:,:,:))./param.mass(2));   
end
%%
Np=[]; Tep=[]; Tip=[]; Mp=[]; Bphip=[]; Bzp=[]; Brp=[]; Rp=[]; Zp=[]; PHIp=[];
nbuf=2;
for izone=1:length(fields)
    a=squeeze(fields(izone).n(1,end,1,1+nbuf:end-nbuf,1+nbuf:end-nbuf));
    Np=[Np a(:)'];
    a=squeeze(mesh(izone).xcells.*param.a0)+param.R0;
    Rp=[Rp a(:)'];
    a=squeeze(mesh(izone).ycells.*param.a0);
    Zp=[Zp a(:)'];
    a=squeeze(fields(izone).T(1,end,1,1+nbuf:end-nbuf,1+nbuf:end-nbuf));
    Tep=[Tep a(:)'];
    a=squeeze(fields(izone).T(2,end,1,1+nbuf:end-nbuf,1+nbuf:end-nbuf));
    Tip=[Tip a(:)'];
    a=squeeze(fields(izone).M(1,end,1,1+nbuf:end-nbuf,1+nbuf:end-nbuf));
    Mp=[Mp a(:)'];
    a=squeeze(Bfield(izone).Bphi);
    Bphip=[Bphip a(:)'];
     a=squeeze(Bfield(izone).Bz);
    Bzp=[Bzp a(:)'];
     a=squeeze(Bfield(izone).Bphi);
    Brp=[Brp a(:)'];

end
Np=Np.*param.n0;
Tep=Tep.*param.T0;
Tip=Tip.*param.T0;
figure;
scatter(Rp,Zp,20,Mp,'filled'); axis equal; colorbar; hold on
plot(mesh(1).wallR,mesh(1).wallZ,'k','linewidth',2)
Bp=sqrt(Brp.^2+Bzp.^2+Bphip.^2);
Mrp=Mp.*Brp./Bp;
Mzp=Mp.*Bzp./Bp;
Mphip=Mp.*Bphip./Bp;
PHIp=ones(size(Rp)).*0.0;
%%

Ztot=size(Np,1);

Nr=250; Nz=Nr*(max(mesh(1).wallZ)-min(mesh(1).wallZ))./(max(mesh(1).wallR)-min(mesh(1).wallR)); Nz=round(Nz);
%Nr=250; Nth=250;
%Nr=99; Nth=101;
clearance=0.01;
R=linspace(min(mesh(1).wallR)-clearance,max(mesh(1).wallR)+clearance,Nr);
Z=linspace(min(mesh(1).wallZ)-clearance,max(mesh(1).wallZ)+clearance,Nz);
hold on; plot([min(R) max(R) max(R) min(R) min(R)],[min(Z) min(Z) max(Z) max(Z) min(Z)],'k')
% phi=linspace(min(PHIp),max(PHIp),32);
% phiu=unique(PHI1); phi=phiu;
% 
% phi=linspace(-20,20,32);% for FT grid
% phiu=unique(PHI1); 
% cond=(phiu>-25).*(phiu<25);
% phiu=phiu(cond>0);
phiu=0;

[R Z]=meshgrid(R,Z);

clear fieldN fieldTe fieldTi fieldM fieldNn fieldTn fieldMI fieldPrad
%% create mask to not extrapolate outside of simulation domain
mask=NaN.*repmat(R,1,1,length(phiu));
for k=1:length(phiu)
    innbouR=[squeeze(mesh(1).xnodes(k,:,1)),squeeze(mesh(2).xnodes(k,:,1))].*param.a0+param.R0;
    innbouZ=[squeeze(mesh(1).ynodes(k,:,1)),squeeze(mesh(2).ynodes(k,:,1))].*param.a0;
   mask(:,:,k)=(1-inpolygon(R,Z,innbouR,innbouZ)).*(inpolygon(R,Z,mesh(1).wallR,mesh(1).wallZ));
   mask(mask==0)=NaN;
end
%%
disp('### Density')

for i=1%:Ztot
    disp(['...density ',num2str(i),'/',num2str(size(Np,1))])
    disp('creating scatteredinterpolant, interpolating...')
    for k=1:length(phiu)
    cond=PHIp==phiu(k);
    n=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Np(i,cond>0)','linear','none') ;
    fieldN(i,1,:,:,k)=n(R,Z)'.*mask(:,:,k)';
    clear n
    end
end
%clear Np
%%
disp('### Mach number')

disp('creating scatteredinterpolant, interpolating...')
    for k=1:length(phiu)
    cond=PHIp==phiu(k);
    mR=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Mrp(1,cond>0)') ;
    mZ=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Mzp(1,cond>0)') ;
    mphi=scatteredInterpolant(Rp(cond>0)',Zp(cond>0)',Mphip(1,cond>0)') ;
    fieldmR(1,1,:,:,k)=mR(R,Z)'.*mask(:,:,k)';
    fieldmZ(1,1,:,:,k)=mZ(R,Z)'.*mask(:,:,k)';
    fieldmphi(1,1,:,:,k)=mphi(R,Z)'.*mask(:,:,k)';
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

% %%
% if exist('Npn')
% disp('### n_N')
% disp('creating scatteredinterpolant...')
% nn=scatteredInterpolant(Rpn',Zpn',PHIpn',Npn') ;
% clear Npn
% disp('interpolating...')
% for iphi=1:length(phi)
%     %disp(['interpolating phi ',num2str(iphi),'/',num2str(length(phi))])
%     fieldNn(1,1,:,:,iphi)=nn(R,Z,ones(size(Z)).*phi(iphi))'.*mask(:,:,iphi)';
%     end
%     clear nn
% %%
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
% end
% %%
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
save(fnameout,'mesh','fieldN','fieldTe','fieldTi','fieldmR','fieldmZ','fieldmphi','fieldMI','fieldNn','fieldTn','-v7.3')
disp(['saved in ',fnameout])



disp('hop')
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
n0=param.n0;
T0=param.T0;
B0=param.R0;
param.R0=1.8;
save(fnameout,'n0','T0','B0','param','-append','-v7.3')

%% save wall
fnamewall= [folder_name{1},'/WEST_wall.txt'];
wall(:,1)=mesh(1).wallR.*100;% in cm
wall(:,2)=mesh(1).wallZ.*100;
save(fnamewall,'wall','-ascii')
disp(['wall saved in ',fnamewall])
%%
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

