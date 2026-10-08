% open source release September 2026
try
    %initialization
    clear all
    test=0;
    if test  mainfolder='/Users/fnespoli/DIS_release/';
                 inputfile='/Users/fnespoli/DIS_release/input_DIS.txt';
    else mainfolder='/p/emc3-eirene/fnespoli/DIS/'; inputfile='input.txt'
    end
    
    addpath(mainfolder)
    tic;


    FID = fopen(inputfile);
    data = textscan(FID,'%s','delimiter','\n');
    a= string(data{:});
    fclose(FID);
    
    rd=str2num(a(2));
    material=data{1}{4};
    x0i=str2num(a(6));
    v0i=str2num(a(8));
    NKU=str2num(a(10));
    phi_or_t=str2num(a(12));
    if length(phi_or_t)>1 
        alpha_0=phi_or_t(2); phi_or_t=phi_or_t(1); 
        disp(['randomized velocity with std of ',num2str(alpha_0),' deg'])
    elseif phi_or_t==4 alpha_0=10;
    end
    fnamepbg=data{1}{14};
    folder1=data{1}{16};
    folder2=data{1}{18};
    fname0=data{1}{20};
    fname0=[fname0,'_',material,'_',num2str(rd),'m'];
    fnamewall=data{1}{22};
    massrate=str2num(a(24));
    sw=str2num(a(26));
    sw_neu=sw(1); sw_imp=sw(2);
    mainion=data{1}{28};
    %% upload plasma background
    disp('upload plasma background...')
    fnamepbg
    load(fnamepbg);
    fieldNe=fieldN(1,:,:,:,:)./n0;
    fieldNi=fieldN./n0;

    if size(fieldNi,1)==1
        disp('No impurities provided in background plasma. Impurities excluded')
    elseif ~sw_imp
        fieldNi=fieldNi(1,:,:,:,:);%to not load impurities
        disp('Impurities excluded')
    end
    clear fieldN
    fieldTe=fieldTe./T0;
    if exist('fieldTi')
        fieldTi=fieldTi./T0;
    else
        fieldTi=fieldTe;
    end
    if exist('fieldNn')
        fieldNn=fieldNn./n0;
        if exist('fieldTn')
            fieldTn=fieldTn./T0;
        end
    end
    if ~exist('fieldNn')
        disp('No neutral density  provided in background plasma. Neutrals excluded')
    elseif ~exist('fieldTn') && exist('fieldNn')
        disp('No neutral Temperature  provided in background plasma. Neutrals temperature set to 1eV')
        fieldTn=ones(size(fieldNn))./T0;
    elseif ~sw_neu
        clear fieldNn fieldTn;% to not load neutrals
        disp('Neutrals excluded')
    end
    mami=1; %m_atom/m_i
    if size(fieldNi,1)==1 %H
        A=[1]; Zi=[1];
    elseif size(fieldNi,1)==2 %He
        A=[4 4]; Zi=[1 2];
    elseif size(fieldNi,1)==7 %H+C
        A=[1 12 12 12 12 12 12];
        Zi=[1 1 2 3 4 5 6];
   elseif size(fieldNi,1)==43 %D+W
        A=[1 184*ones(1,42)];
        Zi=[1 1:42];
    end
    if strcmp(mainion,'H')
        A(1)=1;
    elseif strcmp(mainion,'D')
        A(1)=2;
    else
        disp('main ion not recongnized/implemented')
    end
if strcmp(fnamewall,'NO')
    disp('No wall file specified, reflections with wall ignored')
else
    try
wall=load(fnamewall);
Rwall=wall(:,1)./100; Zwall=wall(:,2)./100; %wall units in cm
clear wall
    catch
     disp('Can not load specified wall file. Reflections with wall ignored')
     clear Rwall Zwall
    end
end
    %% intialization
    %constants
    epsilon0=8.8542e-12;
    mi=A(1).*1.67e-27;%kg, H
    me=9.1094e-31; %kg
    Tw=300;
    qe=1.6022e-19;
    sigma_SB=5.67e-8; %W/m^2/sec/K^4
    h=6.6261e-34; %J*s Planck constant
    kb=1.3807e-23; %J/K Boltzmann constant
    mime=sqrt(mi/me);
    cs0=9.79e3*sqrt(2*T0/A(1));
    rho_s=1.02*sqrt(T0*A(1))/(B0.*1e4); %m
    E_rec=13.6; %eV, ionization potentia energy (H) for recombination on dust surface
    a=rd/rho_s; %dust size in rho-s units
 
    mesh(1).xcells=mesh(1).xcells./rho_s;
    mesh(1).ycells=mesh(1).ycells./rho_s;

    if length(mesh(1).phicells)>1
    mesh(1).phinodes=mesh(1).phicells-mean(diff(mesh(1).phicells))/2;
    mesh(1).phinodes(end+1)=mesh(1).phinodes(end)+diff(mesh(1).phinodes(1:2));
    else
        mesh(1).phinodes=[0 1];
    end
    if max(mesh(1).phicells)>(2.*pi)
        mesh(1).phinodes=mesh(1).phinodes*pi/180;
        mesh(1).phicells=mesh(1).phicells*pi/180;
    end
   
    x0i=[x0i(1)/rho_s x0i(2)/rho_s x0i(3)];% in rho_s
    X0=x0i(1).*cos(x0i(3));  Y0=x0i(1).*sin(x0i(3)); Z0=x0i(2);
    v0i=v0i./cs0; %in c_s
    v00=sqrt(v0i(1).^2+v0i(2).^2+v0i(3).^2);
    
    if exist('t')  
        t_f=(t-t(1)).*1e-3./(rho_s/cs0); %from t_f in ms to normalized units
        clear t
        dt=round(min(diff(t_f))); 
    else
        dt=round(5e-6/(rho_s/cs0));% dt=5e-6 s, normalized
        t_f=0:dt:size(fieldNe,2).*dt;
    end
    % in this code M is Vp/cs0=Vp/cs*(cs/cs0)
    fieldmR=fieldmR.*sqrt(fieldTe); %sqrt*(Te/T0) but Te is already normalized to T0
    fieldmZ=fieldmZ.*sqrt(fieldTe);
    fieldmphi=fieldmphi.*sqrt(fieldTe);
    if exist('fieldER')
        fieldER=fieldER./T0.*rho_s;
        fieldEZ=fieldEZ./T0.*rho_s;
        fieldEphi=fieldEphi./T0.*rho_s;
    end
    
    if strcmp(material,'C')
        Pvap_A=-40181; Pvap_B=14.8; % For vapor pressure [Smirnov PPCF 2007]
        mI=12*1.67e-27;%kg
        rho_d=2267; %graphite kg/m^3
        cp=1709;%DUMBO J/K/kg
        Wf=4.6;% eV for C [Pigarov 2005]2, or S.C. Jain 1952 1952https://doi.org/10.1098/rspa.1952.0116
        hsub=2.97e7; %%DUMBO hsub=2.9e7 J/kg from DUMBO
        epsilon_d=0.75; %Federici PPCF 2003
        b_see=0.4;% For graphite, Smirnov 2007
        boilT=3925;
        meltT=1e6;
    elseif strcmp(material,'B')
        Pvap_A=-28238; Pvap_B=11.98;%For vapor pressure fitted from [American Insitute of Physics handbook 3rd ed. p.4-298]
        mI=10.8*1.67e-27;%kg
        rho_d=2080; % LIQUID kg/m^3
        rho_d=2340; % SOLID kg/m^3
        cp=11.087/10.811e-3;%J/mol/K--> J/K/kg
        Wf=4.45;% CRC handbook p.12-124
        hsub=508.*1e3/10.811e-3; %508 kJ/mol--> J/kg
        hmel=50.2.*1e3/10.811e-3; %50.2 kJ/mol--> J/kg
        epsilon_d=0.27; % L.A. Koval, Combustion, Explosion, and Shock Waves, Vol. 46, No. 2, pp. 178?182, 2010
        epsilon_d=0.8; % total hemispherical emittance, Claude P. Talley Master thesis 1960 "Preparation of elemental boron and measurement of its emittance at incandescent temperatures"
        b_see=0.5; %General relation from Dekker1981 because can't find B specific data
        boilT=4200;
        meltT=2349;
    elseif strcmp(material,'B4C')
        Pvap_A=-31690; Pvap_B=13.42;%For vapor pressure fitted from [H.E. Robson and P.W. Gilles J.Phys.Chem 1964]
        mI=55.255*1.67e-27;%kg
        rho_d=2520; % kg/m^3
        cp=100/55.255e-3;%J/mol/K--> J/K/kg data averaged from https://webbook.nist.gov/cgi/cbook.cgi?ID=C12069328&Mask=2&Type=JANAFS&Table=on#JANAFS
        Wf=(5.15+5.46)/2;% Surface Review and Letters, Vol. 19, No. 4 (2012) 1250040
        hsub=138.*4.184*1e3/55.255e-3; %508 kcal/mol--> J/kg from [H.E. Robson and P.W. Gilles J.Phys.Chem 1964]
        hmel=(1350+2030)/2.*1e3; %50.2 kJ/kg--> J/kg from https://www.azom.com/properties.aspx?ArticleID=75
        epsilon_d=0.93; %Z. Lv Materials 2023 https://doi.org/10.3390/ma16227212
        b_see=0.5; %General relation from Dekker1981 because can't find B specific data
        boilT=3770;
        meltT=3063;
    elseif strcmp(material,'BN')
        Pvap_A=-25710; Pvap_B=13.87;%For vapor pressure fitted from [Lorenz R, Woolcock J (1928) Z Anorg Chemie 176: 283-304]
        mI=24.82*1.67e-27;%kgfi
        rho_d=2100; %hexagonal, 3450 for c-BN
        rho_d=2250; %from Momentive specs
        cp=19.7/24.82e-3;%J/mol/K--> J/K/kg
        cp=45/24.82e-3;%J/mol/K--> J/K/kg  data averaged from https://webbook.nist.gov/cgi/cbook.cgi?ID=C10043115&Units=SI&Mask=2&Type=JANAFS&Plot=on#JANAFS
        Wf=4.5; %electron affinity http://www.ioffe.ru/SVA/NSM/Semicond/BN/bandstr.html
        hmel=81.*1e3/24.82e-3; % kJ/mol--> J/kg CRC
        hsub=hmel; % since it sublimates?
        epsilon_d=0.85; %In the IR range https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5241
        b_see=0.5; %General relation from Dekker1981 because can't find B specific data
        boilT=3246;
        meltT=1e6;%it sublimates at normal pressure, but melts at higher pressure
    elseif strcmp(material,'Li')
        Pvap_A=-8042; Pvap_B=10.02; % For vapor pressure [Smirnov PPCF 2007]
        mI=6.94*1.67e-27;%kg
        rho_d=534 ; % kg/m^3
        cp= 24.860/6.94e-3;%J/mol/K--> J/K/kg;
        cp= 27/6.94e-3;%J/mol/K--> J/K/kg; AVERAGE between melted and solid states
        Wf=2.93;% CRC handbook p.12-124
        hsub=136.*1e3/6.94e-3; % kJ/mol--> J/kg
        hmel=3.0.*1e3/6.94e-3; % kJ/mol--> J/kg
        epsilon_d=0.04; %Rosenberg FED 2009, and A.S. Pricshivichin DOI: 10.21517/0202-3822-2019-42-2-89-95
        b_see=1.3; %After Vignitchouk2014, Schwarz S A 1990 J. Appl. Phys. 68 2382
        boilT=1603;
        meltT=453.65;
    else
        disp([material, ' :material not recognized. Exiting'])
        if ~test exit; end
    end
    
    C_see=load([mainfolder,'YoungDekkerFit_',material,'.txt']);
    m_d=4*pi/3*(a*rho_s).^3.*rho_d;
    
    %normalization/switches
    k=(3/(4*pi*rho_d))^(2/3)*mi*m_d^(-1/3)*n0*rho_s; %ion drag force
    k1=(3/rho_d)^(2/3)*(4*pi/m_d)^(1/3)*mI*rho_s/cs0;
    kg=rho_s/cs0.^2*9.81;
    kE=(3/4/pi/rho_d/m_d.^2).^(1/3).*(T0/cs0).^2.*(4*pi*epsilon0);
    kL=(3/rho_d/m_d.^2*(4*pi).^2).^(1/3).*(T0.*B0.*rho_s/cs0).*(epsilon0);
    
    sw_th=1;%switch for thermionic emission
    sw_see=1;%switch for seconday electron emission
    x=linspace(-8,3,5000);
    HeavisideX=(x>=0); % heaviside(x)
    HeavisidemX=(x<0); % heaviside(-x)
    
    
    %% averaging 3d turbulence    
    fields_list={'fieldNe','fieldNi','fieldTe','fieldTi','fieldER','fieldEZ','fieldEphi','fieldmR','fieldmZ','fieldmphi'};
    if NKU==Inf %2D averaged fields
        disp('averaging 3D fields to 2D...')
        for i=1:length(fields_list)
            if exist(fields_list{i})
            eval([fields_list{i},'=(nanmean(nanmean(',fields_list{i},',2),5));'])
            end
        end
    elseif abs(NKU)>0
        disp('averaging 3D fields to alter Ku number...')
        for i=1:length(fields_list)
            if exist(fields_list{i})
            eval(['f1=',fields_list{i},'=(nanmean(nanmean(',fields_list{i},',2),5));']);
            eval([fields_list{i},'=(',fields_list{i},'-f1)./(NKU+1)+f1;'])
            end
        end
        clear f1 f11
    end
    
    
    T1=toc;
    s = seconds(T1);
    s.Format = 'hh:mm:ss';
    disp(['Time for initialization ',char(s)])
    disp('-------------------------------------------')
    %%
    disp('computing dust dynamics...')
  
        if phi_or_t==1 % change initial toroidal position
            kkktot=length(mesh(1).phicells);
        elseif phi_or_t==2 && size(fieldNe,2)==1
            kkktot=1;
        elseif phi_or_t==2 || phi_or_t==3  % change initial time, or both, randomly
            kkktot=1000; rng('shuffle')
        elseif phi_or_t==4  % random initial velocity direction
            kkktot=500; rng('shuffle')
        elseif phi_or_t==5 % scan R,phi position where plasma is actually there        
            n1=squeeze(nansum(squeeze(fieldNe),2));
            N=10;% resampling in R because too many points
            ilist=[]; jlist=[];
            for i=1:size(n1,1)
                for j=1:size(n1,2)
                if n1(i,j)>0 && any(mesh.xcells(i,1)==mesh.xcells(1:N:end,1))
                    ilist=[ilist i]; jlist=[jlist j];end
                end
            end
          kkktot=length(ilist); clear N n1  i j 
              
        end
     if (NKU==Inf||test) kkktot=1; end
    
    for kkk=1:kkktot
        if (NKU==Inf||test)
            kkk=1;
        else
            
            fname1=dir([folder1,fname0,'*']);
            if isempty(fname1) kstart1=1;
            else
                for l=1:length(fname1)
                    ll = strfind(fname1(l).name,'_');
                    mm = strfind(fname1(l).name,'.mat');
                    nn(l) = str2num(fname1(l).name(ll(end)+1:mm(1)-1));
                end
                klist=1:kkktot;
                kleft=setdiff(klist,nn);
                if isempty(kleft)
                    disp('we have all the files')
                    break;
                else
                    disp(['first missing file detected: ',num2str(kleft(1))])
                    kkk=kleft(1);
                end
            end
        end
        
        fnameout=[folder1,fname0,'_',num2str(kkk,'%03i'),'.mat'];
        fnameoutnew=fnameout;
        save(fnameout,'kkk'); %place holder file
        tic;
        disp(['Trajectory ',num2str(kkk),'/',num2str(kkktot)]);
        
        solvedyn=1; it=1;
         Nfact=2*pi/(mesh(1).phinodes(end)-mesh(1).phinodes(1)); clear x0 v0
        if phi_or_t==1 %change the toroidal position
            itf0=1; x0=[x0i(1), x0i(2), x0i(3)+(kkk-1)*mean(diff(mesh(1).phicells))]; 
        elseif phi_or_t==2 %change initial time
            if length(t_f)==1 tf0=1; else 
            tf0=rand(1)*(t_f(end)-t_f(1))+t_f(1);
            itf0=dsearchn(t_f',tf0);
            end
        elseif phi_or_t==3 %both
            tf0=rand(1)*(t_f(end)-t_f(1))+t_f(1);
            itf0=dsearchn(t_f',tf0);
            x0=[x0i(1), x0i(2), mesh(1).phicells( round( rand(1)*(length(mesh(1).phicells)-1)+1))];
            disp(['itf0=',num2str(itf0),', iphi0=',num2str(dsearchn(mesh(1).phicells',mod(x0(3),2*pi/Nfact)))])
        elseif phi_or_t==4 %randomize initial velocity direction
            itf0=1; 
            alpha=180+alpha_0.*randn(1);% normally distributed
            beta=rand(1)*2*180-180; %uniformily distributed
            v0(1,1)=v00.*sind(alpha).*cosd(beta);
            v0(2,1)=v00.*cosd(alpha);
            v0(3,1)=v00.*sind(alpha).*sind(beta);       
            X00=X0+5e-3./rho_s.*randn(1);
            Y00=Y0+5e-3./rho_s.*randn(1);
            Z00=Z0;
            if v0i(1)~=0 %case of horizontal injection
                X00=X0;
                Y00=Y0+5e-3./rho_s.*randn(1);
                Z00=Z0+5e-3./rho_s.*randn(1);
                alpha=180+alpha_0.*randn(1);% normally distributed
                v0(1,1)=v00.*cosd(alpha); 
                v0(2,1)=v00.*sind(alpha).*cosd(beta);
                v0(3,1)=v00.*sind(alpha).*sind(beta); 
            if v0i(2)~=0 || v0i(3)~=0 %case of arbitrary direction
                aa=[-1 0 0 ]' ; 
                bb=v0i'./sqrt(v0i(1).^2+v0i(2).^2+v0i(3).^2);
                R=2.*((aa+bb)*(aa+bb)')/((aa+bb)'*(aa+bb))-eye(3);
                v0=R*v0;
                X00=X0; Y00=Y0; Z00=Z0;
            end
            end
            x0=[sqrt(X00.^2+Y00.^2) Z00 atan2(Y00,X00)];
            disp(['x0 = [',num2str(x0),'], v0 = [',num2str(v0'),']'])
        elseif phi_or_t==5 %scan radial and toroidal
            itf0=1;
            x0=[mesh.xcells(ilist(kkk),1), x0i(2), mesh.phicells(jlist(kkk))];
            disp(['i_R = ',num2str(ilist(kkk)),'/',num2str(length(mesh.xcells(:,1))),' i_phi = ',num2str(jlist(kkk)),'/',num2str(length(mesh.phicells))])
        end
        clear t R vR Z vZ phi vphi Chi Tau Td Mu ITHd IRd IZd IPHId Nd Ted Md ERd EZd Ephid Fi FE  dead MTOR LNLAMBDA
        %%
        Ntmax=4e5;
        R=zeros( Ntmax,1).*NaN;
        t=R;
        vR=R;
        Z=R;
        vZ=R;
        phi=R;
        vphi=R;
        Chi=R;
        Tau=zeros( Ntmax,size(fieldNi,1)).*NaN;
        Td=R; Mu=R;
        IPHId=R;
        IRd=R;
        ITHd=R;
        IZd=R;
        Nd=zeros( Ntmax,size(fieldNi,1)+1).*NaN;
        Ted=R;
        Md=zeros( Ntmax,3).*NaN;
        Ed=zeros( Ntmax,3).*NaN;
        Psimelt=R;
        Nnd=R; Tnd=R;         Rho_d = R;  Cp = R;
        %%
        it=1;
        
        while solvedyn
            if it==1
                if ~exist("x0") x0=x0i; end
                if ~exist("v0") v0=v0i; end
                itf=itf0; R(it)=x0(1); vR(it)=v0(1); Z(it)=x0(2); vZ(it)=v0(2); phi(it)=x0(3); vphi(it)=v0(3);
                Chi(it)=3; Tau(it,:)=ones(size(A)); Td(it)=[300]; mu=1;
                Mu(it)=mu;
                IPHId(it)=dsearchn(mesh(1).phicells',x0(3));
                IRd(it)=dsearchn(mesh(1).xcells(:,1),x0(1));
                ITHd(it)=dsearchn(mesh(1).ycells(1,:)',x0(2));
                IZd(it)=1;
                Nd(it,:)=NaN.*ones(size(fieldNi,1)+1,1); Ted(it)=NaN;
                Nnd(it)=NaN; Tnd(it)=NaN;
                Md=[0 0 0]; Ed=[0 0 0];
                
                Fi(it)=NaN; FE(it)=NaN;
                Psimelt(it)=0; psimelt=0;
                t(1)=t_f(itf);
            end
            
            ir=IRd(it); ith=ITHd(it); iz=IZd(it);
            if min(mesh(1).phicells)+max(mesh(1).phicells)==0
                iphi=dsearchn(mesh(1).phicells',mod(phi(it)-pi/Nfact,2*pi/Nfact)-pi./Nfact);
            elseif min(mesh(1).phicells)>=0
                iphi=dsearchn(mesh(1).phicells',mod(phi(it),2*pi/Nfact));
            end
            if NKU~=Inf
                iphi1=iphi;
                itf=dsearchn(t_f',mod(t(it),t_f(end)));
            else
                itf=1; iphi1=1;
            end
            
            ITF(it)=itf;
            Ni=squeeze(fieldNi(:,itf,ir,ith,iphi1))';
            Ne=fieldNe(1,itf,ir,ith,iphi1);
            if Ne<0  dead=-1; break;end
            Te=fieldTe(1,itf,ir,ith,iphi1);
            Ti=squeeze(fieldTi(:,itf,ir,ith,iphi1))';
            if exist('fieldNn')
                Nn=fieldNn(1,itf,ir,ith,iphi1);
                Tn=fieldTn(1,itf,ir,ith,iphi1);
                Mn=0;
            else
                Nn=0; Tn=0; Mn=0;
            end
           
       
            if exist('fieldER')
                ER=fieldER(1,itf,ir,ith,iphi1);
                EZ=fieldEZ(1,itf,ir,ith,iphi1);
                Ephi=fieldEphi(1,itf,ir,ith,iphi1);
            else
                ER=0; EZ=0; Ephi=0;
            end
            
            mR=fieldmR(iz,itf,ir,ith,iphi1);
            mZ=fieldmZ(iz,itf,ir,ith,iphi1);
            mphi=fieldmphi(iz,itf,ir,ith,iphi1);
            %temperature dependent density and thermal capacity
           [rho_d cp]=cp_func(material, Td(it), meltT);
                if (isnan(Ne) || Te.*T0<0.1)  mR=0; mZ=0; mphi=0;  Nn=0; Tn=0; Ne=1e13/n0; Te=0.1/T0; Ti=Te; Ni=Ne.*ones(1,length(Zi)); 
                end

                if length(Ti)==1 Ti=Ti.*ones(size(Ni)); end
                tau=Ti./Te;
                u=sqrt((mR-vR(it)).^2+(mZ-vZ(it)).^2+(mphi-vphi(it)).^2)./sqrt(Ti).*sqrt(A./A(1));
                lD=debyelength(Te.*T0,Ne.*n0);
                ad=(3/4/pi*m_d*mu/rho_d).^(1/3);
                lnLambda=Coul_log(Chi(it),tau,ad,u,lD);
                %% dust charge
                % all currents are divided by pi*a_d^2*e*cs0*n_e
                
                T=Td(it);
                %electrons
                Iep=-sqrt(Te.*4/pi).*mime.*exp(-x);
                Ien=-sqrt(Te.*4/pi).*mime.*(1-x);
                %ions
                Iip=0; Iin=0;
                for z=1:length(Zi)
                    Iip=Iip + Ni(z)./Ne.*Zi(z).*sqrt(Ti(z).*A(1)./A(z))*(erf(u(z)).*(u(z)+1/(2*u(z))+x.*Zi(z)/(tau(z)*u(z)))+1/sqrt(pi)*exp(-u(z).^2));
                    um=u(z)-sqrt(Zi(z).*abs(x)./tau(z));
                    up=u(z)+sqrt(Zi(z).*abs(x)./tau(z));
                    Iin=Iin + Ni(z)./Ne.*Zi(z).*sqrt(Ti(z).*A(1)./A(z))/4/u(z).*((erf(up(z))+erf(um(z))).*(1+2*u(z).^2+2*x.*Zi(z)/tau(z))+2/sqrt(pi)*(exp(-(um(z)).^2).*(up(z))+exp(-(up(z)).^2).*(um(z))));
                end
                %thermionic emission
                Ithp=sw_th.*16.*pi.*me*(kb.*T).^2./h.^3./cs0./Ne./n0.*exp(-Wf.*qe./kb./T);
                Ithn=Ithp.*(1-x.*(Te.*T0.*qe)./(T.*kb)).*exp(x.*(Te.*T0.*qe)./(T.*kb));
                Ithn(x>0)=0;
                %SEE emission
                F_f = 10^sum(C_see.*log10(Te.*T0).^((0:numel(C_see)-1)'));
                dSEEp=sw_see.*2/(2-b_see)*F_f;
                dSEEn=sw_see.*2/(2-b_see)*F_f*(1-3*x.*Te.*T0./Wf)./(1-x)./(1-x.*Te.*T0./Wf).^3.*exp(-x);%Autricque These
                
                eqtot=(Iep.*(1-dSEEp)+Iip+Ithp).*HeavisideX+(Ien.*(1-dSEEn)+Iin+Ithn).*HeavisidemX;
                
                %% dust temperature
                %all fluxes are divided by 4pia^2
                
                %chi>=0
                q_currp=Ne.*Te.*abs(Iep).*(2+x);  %Autrique these, electron contribution
                %chi<0
                q_currn=Ne.*Te.*abs(Ien)./(1-x).*(2-x);
                for z=1:length(A) %sum contribution for various ions
                    q_currp=q_currp+ ...
                        Ni(z).* sqrt(Ti(z).*A(1)./A(z)).*Ti(z).*1/4.*(2/sqrt(pi).*(5+2.*u(z).^2+2.*x.*Zi(z)./tau(z)).*exp(-u(z).^2)... %Smirnov PPCF 2007
                        +(3./u(z)+12.*u(z)+4.*u(z).^3+2.*x.*Zi(z)/tau(z).*(1./u(z)+2.*u(z))).*erf(u(z)));
                    um=u(z)-sqrt(Zi(z).*abs(x)./tau(z));
                    up=u(z)+sqrt(Zi(z).*abs(x)./tau(z));
                    q_currn=q_currn + ...
                        Ni(z).* sqrt(Ti(z).*A(1)./A(z)).*Ti(z).*1/8.*(2/sqrt(pi).*((5+2.*u(z).^2-(3+2.*u(z).^2)./u(z).*sqrt(-x.*Zi(z)/tau(z))).*exp(-up(z).^2)...
                        +(5+2.*u(z).^2+(3+2.*u(z).^2)./u(z).*sqrt(-x.*Zi(z)/tau(z))).*exp(-um(z).^2))...
                        +(3./u(z)+12.*u(z)+4.*u(z).^3+2.*x.*Zi(z)/tau(z).*(1./u(z)+2.*u(z))).*(erf(up(z))+erf(um(z))));
                end
                %energy of recombination on dust surface
                %This is now holding just in the case of hydrogen, no impurities
                q_currp= q_currp+E_rec/T0.*abs(Iip).*Ni(1);
                q_currn= q_currn+E_rec/T0.*abs(Iin).*Ni(1);
                q_currp=(q_currp).*n0.*cs0.*T0./4;
                q_currn=(q_currn).*n0.*cs0.*T0./4;
                
                q_seep=3.*Wf./T0.*Ne.*Te.*abs(dSEEp.*Iep) .*n0.*cs0.*T0./4;%Smirnov PPCF 2007
                zeta=-x.*Te.*T0./Wf;
                q_seen=3.*Wf./T0.*Ne.*Te.*abs(dSEEn.*Ien) .*n0.*cs0.*T0./4 .*(1+zeta).*(1+2.*zeta)./(1+3.*zeta);%Vignitchouk PPCF 2014
                
                
                qred=(q_currp-q_seep).*HeavisideX+(q_currn-q_seen).*HeavisidemX; %heat flux/4pia^2 excluding the terms depending on Td

                %% solve for charge and temperature 

                i=find(diff(sign(eqtot)));
                x1=(x(1:end-1)+x(2:end))/2;
                chi=x1(i);
                if length(chi)>1 
                    j=dsearchn(chi',Chi(it));
                    i=i(j); chi=chi(j);
                end 
                qred=interp1(x,qred,chi);
                
                ktemp=(4.*pi/mu/m_d)^(1/3).*(3/rho_d).^(2/3)./cp.*(rho_s./cs0).*qe;
                [ts,td] = ode45(@(t,y) eqTemp(t,y,ktemp,qred,chi,sw_th,Te.*T0,Wf,Pvap_A,Pvap_B,hsub,mI,epsilon_d,Tw,kb,qe,me,h,sigma_SB),[0,dt/2,dt],[Td(it)]);
                td=td(end);
                if ~isreal(td) 
                disp('Dust temperature is not real')
                if test keyboard; else dead=-10; break; end
                end
                if Td(it)<=meltT && td>=meltT
                    disp('melting')
                    qmelt=(m_d*mu/4/pi).^(2/3)*(rho_d/3)^1/3*hmel*(1-Psimelt(it));
                    [ts,td] = ode45(@(t,y) eqTemp(t,y,ktemp,qred-qmelt,chi,sw_th,Te.*T0,Wf,Pvap_A,Pvap_B,hsub,mI,epsilon_d,Tw,kb,qe,me,h,sigma_SB),[0,dt/2,dt],[Td(it)]);
                    td=td(end); psimelt=1;
                    if td<=meltT
                        td=meltT;
                        q0=qred-q_sub(Td(it),Pvap_A,Pvap_B,hsub,mI,kb,qe)-sw_th.*q_th(Td(it),chi,Wf,Te,me,kb,h,qe)-q_bb(Td(it),Tw,epsilon_d,sigma_SB);
                        q1=qred-q_sub(td,Pvap_A,Pvap_B,hsub,mI,kb,qe)-sw_th.*q_th(td,chi,Wf,Te,me,kb,h,qe)-q_bb(Td,Tw,epsilon_d,sigma_SB);
                        q2=q0-q1;
                        psimelt=Psimelt(it)+q2.*(1/hmel)*(3/rho_d).^(2/3).*(4*pi/mu/m_d).^(1/3).*dt.*(rho_s/cs0);
                    end
                elseif Td(it)>=meltT && td<=meltT
                    disp('solidification')
                    qmelt=(m_d*mu/4/pi).^(2/3)*(rho_d/3)^1/3*hmel*(Psimelt(it));
                    [ts,td] = ode45(@(t,y) eqTemp(t,y,ktemp,qred-qmelt,chi,sw_th,Te.*T0,Wf,Pvap_A,Pvap_B,hsub,mI,epsilon_d,Tw,kb,qe,me,h,sigma_SB),[0,dt/2,dt],[Td(it)]);
                    td=td(end); psimelt=0;
                    if td>=meltT
                        td=meltT;
                        q0=qred-q_sub(Td(it),Pvap_A,Pvap_B,hsub,mI,kb,qe)-sw_th.*q_th(Td(it),chi,Wf,Te,me,kb,h,qe)-q_bb(Td(it),Tw,epsilon_d,sigma_SB);
                        q1=qred-q_sub(td,Pvap_A,Pvap_B,hsub,mI,kb,qe)-sw_th.*q_th(td,chi,Wf,Te,me,kb,h,qe)-q_bb(Td,Tw,epsilon_d,sigma_SB);
                        q2=q0-q1;
                        psimelt=Psimelt(it)+q2.*(1/hmel)*(3/rho_d).^(2/3).*(4*pi/mu/m_d).^(1/3).*dt.*(rho_s/cs0);
                    end
                end
                
                %% dust mass
                [ts,mu] = ode45(@(t,y) eqMass(t,y,td,k1,Pvap_A,Pvap_B,mI,kb),[0,dt/2,dt],[mu]);
                mu=mu(end);
                if mu<=1e-4 disp('dust evaporated completely'); dead=4; break; end
                
                
           % end
            %% dust trajectory
            %if it==17558 keyboard; end
            [ts,Rs] = ode45(@(t,y) eqmotR(t,y,Ni,mR,k,mu,Ti,vphi(it),chi,lnLambda,Te,ER,kE,A,Zi,Nn,Tn,Mn,mami),[0,dt/2,dt],[R(it) vR(it)]);
            [ts,Zs] = ode45(@(t,y) eqmotZ(t,y,Ni,mZ,k,mu,Ti,chi,lnLambda,kg,Te,EZ,kE,A,Zi,Nn,Tn,Mn,mami),[0,dt/2,dt],[Z(it) vZ(it)]);
            [ts,phis] = ode45(@(t,y) eqmotphi(t,y,Ni,mphi,k,mu,Ti,vphi(it),vR(it),R(it),chi,lnLambda,Te,Ephi,kE,A,Zi,Nn,Tn,Mn,mami),[0,dt/2,dt],[phi(it) vphi(it)]);
            
            it=it+1; 
            t(it)=t(it-1)+dt;
            R(it)= Rs(end,1);
            vR(it)= Rs(end,2);
            Z(it)= Zs(end,1);
            vZ(it)= Zs(end,2);
            phi(it)= phis(end,1);
            vphi(it)=phis(end,2);
            Chi(it)= chi;
            Tau(it,:)= tau;
            Td(it)= td; Mu(it)= mu;
            IPHId(it)=dsearchn(mesh(1).phicells',phi(it));
            IRd(it)=dsearchn(mesh(1).xcells(:,1),R(it));
            ITHd(it)=dsearchn(mesh(1).ycells(1,:)',Z(it));
            IZd(it)=1;
            Nd(it,:)= [Ne Ni]; Ted(it)= Te;
            Md(it,:)= [mR mZ mphi];
            Ed(it,:)= [ER EZ Ephi];
            Psimelt(it)=psimelt;
            Nnd(it)=Nn; Tnd(it)=Tn;
            Rho_d(it) = rho_d; Cp(it) = cp;
            
            if exist('Rwall')
            [R,Z,vR,vZ]=wall_reflection(Rwall./rho_s,Zwall./rho_s,R,Z,vR,vZ,it);
            end
            
            if IZd(it)==2 && (ITHd(it)-ITHd(it-1))<(-mesh(1).Ntheta/4)
                disp('dust collided with the limiter')
                dead=0;
                break
            end
            if it>Ntmax dead=-1; break;end
            if ~exist('Rwall')
            if R(it)>max(mesh(1).xcells(:,1)) || R(it)<min(mesh(1).xcells(:,1)) || Z(it)>max(mesh(1).ycells(1,:)) || Z(it)<min(mesh(1).ycells(1,:))
                disp('dust outside simulation domain')
                dead=1;
                break
            end
            end
        end
        
        t=t(1:it);
        R=R(1:it);
        vR=vR(1:it);
        Z=Z(1:it);
        vZ=vZ(1:it);
        phi=phi(1:it);
        vphi=vphi(1:it);
        Chi=Chi(1:it);
        Tau=Tau(1:it,:);
        Td=Td(1:it); Mu=Mu(1:it);
        IPHId=IPHId(1:it);
        IRd=IRd(1:it);
        ITHd=ITHd(1:it);
        IZd=IZd(1:it);
        Nd=Nd(1:it,:);
        Ted=Ted(1:it);
        Md=Md(1:it,:);
        Ed=Ed(1:it,:);
        Psimelt=Psimelt(1:it);
        Nnd=Nnd(1:it); Tnd=Tnd(1:it);
        Rho_d=Rho_d(1:it); Cp=Cp(1:it);
        fnameout=[folder1,fname0,'_',num2str(kkk,'%03i'),'.mat'];
        
        save(fnameout,'NKU','n0','T0','B0','cs0','rho_s','a','x0','v0','R','vR','Z','vZ','phi','vphi','Chi','Tau','Td','Mu',...
            'ITHd','IRd','IZd','IPHId','Nd','Ted','Md','Ed','material','t','dt','test','dead','Psimelt','Nnd','Tnd','Rho_d','Cp');
        disp(['results saved in ',fnameout])
        clear fnameoutnew
        T1(kkk+1)=toc;
        s = seconds(T1(kkk+1));
        s.Format = 'hh:mm:ss';
        disp(['Computation time for this step: ',char(s)])
        disp('-------------------------------------------')
    end
    tic;
    fname1=dir([folder1,fname0,'*']);
    if length(fname1)==kkktot         
        if ~exist('fnameout_f')
            fnameout_f=[folder1,'ongoing_',fname0,'.txt'];
            save(fnameout_f,'fnameout_f')% place holder
        reassamble_dust_results(folder1,folder2,fname0)
        disp('Computing neutral atom source...')
        load([folder2,fname0,'.mat'])
        NAS=NeutralAtomSource(rd,Rho_d,mI/mi,t.*rho_s/cs0,Mu);
        save([folder2,fname0,'.mat'],'NAS','-append')
        fnameBsource=[folder2,fname0,'_NAS_EMC3.txt'];
        a=NAS_to_EMC3(R*rho_s*100,Z*rho_s*100,rad2deg(phi),NAS,massrate,mI/(mi./A(1)),fnameBsource,Mu);
        if strcmp(material,'B') || strcmp(material,'C') || strcmp(material,'Li')
            fname_imp=[folder2,fname0,'_input.imp_',material];
        write_EMC3_imp_input(a,fname_imp,material)
        elseif strcmp(material,'B4C') || strcmp(material,'BN')
            disp('rescaling atom source for composite material')
            [m1 m2 a b]=rescale_NAS_composites(fnameBsource,material);
            fname_imp1=[folder2,fname0,'_input.imp_',m1];
            fname_imp2=[folder2,fname0,'_input.imp_',m2];
            write_EMC3_imp_input(a,fname_imp1,m1)
            write_EMC3_imp_input(b,fname_imp2,m2)
        end
        eval(['!rm ',fnameout_f]);
        else
            disp('reassembling of trajectory already ongoing ')
        end
    else
        disp('total number of file not enough for reassembling')
    end
    
    T1(end+1)=toc;
    s = seconds(sum(T1));
    s.Format = 'hh:mm:ss';
    disp(['total computation time ',char(s)])
    if ~test exit; end
    
catch e
    try e
        e.stack
    catch; end
    disp(e.identifier)
    disp(e.message)
    disp(['file', e.stack.file, ' line ', num2str(e.stack.line)])
    disp('exiting')
    if exist('fnameoutnew') eval(['!rm ',fnameoutnew]); end 
    if ~test exit; end
end

%% visualization
if test
    kphi=dsearchn(mesh(1).phicells',x0(3)); dph=0;
    figure; pcolor(mesh(1).xcells.*rho_s,mesh(1).ycells.*rho_s,(squeeze(fieldmphi(1,1,:,:,kphi)))); shading flat; axis equal
%figure; pcolor(mesh(1).xcells.*rho_s,mesh(1).ycells.*rho_s,(squeeze(sum(fieldNi(2:end,1,:,:,kphi),1))).*n0); shading flat; axis equal

    hold on; plot(R'.*rho_s,Z'.*rho_s,'r','linewidth',2)
                if exist('Rwall')
                    plot(Rwall,Zwall,'k','linewidth',2)
                end
end



