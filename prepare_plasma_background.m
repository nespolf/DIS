
try
%initialization
clear all
mainfolder='/p/emc3-eirene/fnespoli/'
%    mainfolder='/Users/fnespoli/'
    addpath ([mainfolder,'scripts'])
    addpath ([mainfolder,'scripts/dust_dyn'])
tic;  


curdir=pwd;
cd ([mainfolder,'TK3X'])
for i=1:10
    folder_list{i}= ['TCVlim_',num2str(20+i),'/Output/'];
end

param = tk3x_load_param(folder_list{1});
mesh = tk3x_load_mesh(folder_list{1},param,0);
Bfield = tk3x_load_Bfield(folder_list{1},mesh,param.R0,0);


IR2=[]; IZ2=[]; ITH2=[];  X=[]; Y=[];
for iz=1:length(mesh)
    [IR1 ITH1]=meshgrid(1:mesh(iz).Nr,1:mesh(iz).Ntheta);
    IR1=double(IR1)'; ITH1=double(ITH1)'; IZ1=ones(size(IR1)).*double(iz);
    if iz>1 IR1=IR1+double(mesh(mesh(iz).neighb(1)).Nr); end
    IR2=[IR2 reshape(IR1,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
    ITH2=[ITH2 reshape(ITH1,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
    IZ2=[IZ2 reshape(IZ1,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
    X=[X reshape(mesh(iz).xcells,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
    Y=[Y reshape(mesh(iz).ycells,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
end
IR=scatteredInterpolant(X',Y',IR2','linear') ;
ITH=scatteredInterpolant(X',Y',ITH2','linear') ;
IZ=scatteredInterpolant(X',Y',IZ2','linear') ;
clear IR2 ITH2 IZ2 ITH1 IR1 IZ1
ator2=[]; aR2=[]; aZ2=[];
for iz=1:length(mesh)
    Bfield(iz).B=sqrt(Bfield(iz).Btor.^2+Bfield(iz).Bpol.^2);
    dR=(mesh(iz).xnodes(1:end-1,2:end)+mesh(iz).xnodes(2:end,2:end))/2-(mesh(iz).xnodes(1:end-1,1:end-1)+mesh(iz).xnodes(2:end,1:end-1))/2;
    dZ=(mesh(iz).ynodes(1:end-1,2:end)+mesh(iz).ynodes(2:end,2:end))/2-(mesh(iz).ynodes(1:end-1,1:end-1)+mesh(iz).ynodes(2:end,1:end-1))/2;
    Bfield(iz).BR=Bfield(iz).Bpol.*dR./sqrt(dR.^2+dZ.^2);
    Bfield(iz).BZ=Bfield(iz).Bpol.*dZ./sqrt(dR.^2+dZ.^2);
    ator1=(Bfield(iz).Btor./Bfield(iz).B);
    aR1=(Bfield(iz).BR./Bfield(iz).B);
    aZ1=(Bfield(iz).BZ./Bfield(iz).B);
    ator2=[ator2 reshape(ator1,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
    aR2=[aR2 reshape(aR1,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
    aZ2=[aZ2 reshape(aZ1,1,mesh(iz).Nr.*mesh(iz).Ntheta)];
end
ator=scatteredInterpolant(X',Y',ator2','linear') ;
aR=scatteredInterpolant(X',Y',aR2','linear') ;
aZ=scatteredInterpolant(X',Y',aZ2','linear') ;
clear aR1 aZ1 ator1 ator2 aR2 aZ2 X Y

%% upload plasma background

 disp('loading fields...')

    fields_list = {'N','Te','Ti','Gammai','PHI'};%,'W'}%,'Ti','Te'}%,'Gammai','Gammae'}%,'Gammai'}%,'Te','Ti'}%,'PHI','Gammai','Gammae'};%,'W','Gammai','Te','Ti'};
fields_list = {'N','Te','Ti'};%
    Ntload_force = -1;plotbool =0; dtmin = 0;

            Nt0=0;
            dtload=4;
            [~,~,~,~,tt,fields] = tk3x_load_all(folder_list,fields_list,Nt0,dtload,plotbool,Ntload_force);

%%

itlast=0;

T1=toc;
disp([num2str(T1),' seconds for inizialization'])
folder1=[mainfolder,'temp/'];
fname0='TK3X_3Dt_plasma_background';       

        
            fname1=dir([folder1,fname0,'*']);
            if isempty(fname1) itstart=1;
            else itstart=length(fname1)+1;
                for l=1:length(fname1) dtnm(l)=fname1(l).datenum; end
                [lll lm]=max(dtnm);
                disp(['last file detected ',fname1(lm).name])
                disp(['restarting from k=',num2str(itstart)])
            end
      
 disp('computing E, EXB...')
 if itstart<size(fields(1).N,1)
for itf=itstart:size(fields(1).N,1)
    tic;
    disp(['it=',num2str(itf),'/',num2str(length(tt))])
    ER=zeros(length(mesh),mesh(1).Nr,mesh(1).Ntheta,mesh(1).Nphi).*NaN;
    EZ=ER; Ephi=ER; mR=ER; mZ=ER; mphi=ER;
 for iz=1:length(mesh)
     for ir=1:mesh(iz).Nr
         for ith=1:mesh(iz).Ntheta
             for iphi=1:mesh(iz).Nphi

    if ir==1
        if mesh(iz).neighb(1)<1
            Er=0;
        else
            Er=-(fields(iz).PHI(itf,ir+1,ith,iphi)-fields(mesh(iz).neighb(1)).PHI(itf,end,ith,iphi))./sqrt((mesh(iz).xcells(ir+1,ith)-mesh(mesh(iz).neighb(1)).xcells(end,ith)).^2+(mesh(iz).ycells(ir+1,ith)-mesh(mesh(iz).neighb(1)).ycells(end,ith)).^2);
        end
    elseif ir==mesh(iz).Nr
        if mesh(iz).neighb(2)<1
            Er=0;
        else
            Er=-(fields(mesh(iz).neighb(2)).PHI(itf,1,ith,iphi)-fields(iz).PHI(itf,ir-1,ith,iphi))./sqrt((mesh(mesh(iz).neighb(2)).xcells(1,ith)-mesh(iz).xcells(ir-1,ith)).^2+(mesh(mesh(iz).neighb(2)).ycells(1,ith)-mesh(iz).ycells(ir-1,ith)).^2);
        end
    else
        Er=-(fields(iz).PHI(itf,ir+1,ith,iphi)-fields(iz).PHI(itf,ir-1,ith,iphi))./sqrt((mesh(iz).xcells(ir+1,ith)-mesh(iz).xcells(ir-1,ith)).^2+(mesh(iz).ycells(ir+1,ith)-mesh(iz).ycells(ir-1,ith)).^2);
    end
    if ith==1 || ith==mesh(iz).Ntheta
        Eth=0;
    else
        Eth=-(fields(iz).PHI(itf,ir,ith+1,iphi)-fields(iz).PHI(itf,ir,ith-1,iphi))./sqrt((mesh(iz).xcells(ir,ith+1)-mesh(iz).xcells(ir,ith-1)).^2+(mesh(iz).ycells(ir,ith+1)-mesh(iz).ycells(ir,ith-1)).^2);
    end
    
        if iphi>1 && iphi<mesh(iz).Nphi
            Ephi1=-(fields(iz).PHI(itf,ir,ith,iphi+1)-fields(iz).PHI(itf,ir,ith,iphi-1))./(mesh(iz).xcells(ir,ith).*2.*mean(diff(mesh(iz).phinodes)));
        elseif iphi==1
            Ephi1=-(fields(iz).PHI(itf,ir,ith,iphi+1)-fields(iz).PHI(itf,ir,ith,end))./(mesh(iz).xcells(ir,ith).*2.*mean(diff(mesh(iz).phinodes)));
        elseif iphi==mesh(iz).Nphi
            Ephi1=-(fields(iz).PHI(itf,ir,ith,1)-fields(iz).PHI(itf,ir,ith,end))./(mesh(iz).xcells(ir,ith).*2.*mean(diff(mesh(iz).phinodes)));
        end
        
    Ephi(iz,ir,ith,iphi)=Ephi1;
    theta=atan2(mesh(iz).ycells(ir,ith),mesh(iz).xcells(ir,ith)-param.R0);
    ER(iz,ir,ith,iphi)=Er.*cos(theta)-Eth*sin(theta);
    EZ(iz,ir,ith,iphi)=Er.*sin(theta)+Eth*sin(theta+pi/2);
    
    EXBR=-(EZ(iz,ir,ith,iphi)*Bfield(iz).Btor(ir,ith)-Ephi(iz,ir,ith,iphi)*Bfield(iz).BZ(ir,ith))./Bfield(iz).B(ir,ith).^2;
    EXBZ=-(-ER(iz,ir,ith,iphi)*Bfield(iz).Btor(ir,ith)+Ephi(iz,ir,ith,iphi)*Bfield(iz).BR(ir,ith))./Bfield(iz).B(ir,ith).^2;
    EXBphi=-(ER(iz,ir,ith,iphi)*Bfield(iz).BZ(ir,ith)-EZ(iz,ir,ith,iphi)*Bfield(iz).BR(ir,ith))./Bfield(iz).B(ir,ith).^2;

    R=mesh(iz).xcells(ir,ith);
    Z=mesh(iz).ycells(ir,ith);
    M=fields(iz).Mi(itf,ir,ith,iphi); 
    
    mphi(iz,ir,ith,iphi)=M.*(ator(R,Z))+EXBphi;
    mR(iz,ir,ith,iphi)=M.*(aR(R,Z))+EXBR;
    mZ(iz,ir,ith,iphi)=M.*(aZ(R,Z))+EXBZ;
    
%    fields(iz).ER(itf,ir,ith,iphi)=ER;
%    fields(iz).EZ(itf,ir,ith,iphi)=EZ;
%    fields(iz).Ephi(itf,ir,ith,iphi)=Ephi;
%    fields(iz).MR(itf,ir,ith,iphi)=mR;
%    fields(iz).MZ(itf,ir,ith,iphi)=mZ;
%    fields(iz).Mphi(itf,ir,ith,iphi)=mtor;
%    itlast=itf;
             end
         end
     end
 end
 T1(end+1)=toc;
 disp('saving...')
 fname=[folder1,fname0,'_',num2str(itf),'.mat']
save(fname,'ER','EZ','Ephi','mR','mZ','mphi','-v7.3')

disp(['computation time ',num2str(T1(end)),' s'])
end
 end

disp('reassembling...')

    fieldER=zeros(length(mesh),size(fields(1).N,1),mesh(1).Nr,mesh(1).Ntheta,mesh(1).Nphi).*NaN;
    fieldEZ=fieldER; fieldEphi=fieldER; fieldmR=fieldER; fieldmZ=fieldER; fieldmphi=fieldER;
    fieldN=fieldER; fieldTe=fieldER; fieldTi=fieldER;
for iz=1:length(mesh)   
    fieldN(iz,:,:,:,:)=fields(iz).N;
    fieldTe(iz,:,:,:,:)=fields(iz).Te;
    fieldTi(iz,:,:,:,:)=fields(iz).Ti;
end
clear fields
for itf=1:size(fieldN,2)
  fname=[folder1,fname0,'_',num2str(itf),'.mat']
load(fname)
for iz=1:length(mesh)
    fieldER(iz,itf,:,:,:)=squeeze(ER(iz,:,:,:));
    fieldEZ(iz,itf,:,:,:)=squeeze(EZ(iz,:,:,:));
    fieldEphi(iz,itf,:,:,:)=squeeze(Ephi(iz,:,:,:));
    fieldmR(iz,itf,:,:,:)=squeeze(mR(iz,:,:,:));
    fieldmZ(iz,itf,:,:,:)=squeeze(mZ(iz,:,:,:));
    fieldmphi(iz,itf,:,:,:)=squeeze(mphi(iz,:,:,:));
end

end


fnamesave=([mainfolder,'TK3X/',fname0,'.mat'])
save(fnamesave, 'fieldN','fieldTe','fieldTi','fieldER','fieldEZ','fieldEphi','fieldmR','fieldmZ','fieldmphi','-v7.3')
T1(end+1)=toc;
disp(['total computation time ',num2str(sum(T1)),' s'])
exit

catch e
    try disp(e.stack); catch; end
    disp(e.identifier)
    disp(e.message)
    disp('exiting')
    exit
end
