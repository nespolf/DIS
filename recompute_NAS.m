try 
clear all
    test=1;
    if test  mainfolder='/Users/fnespoli/scripts/DIS/';
         %      inputfile='/Users/fnespoli/run_dis/input_WEST.txt'
%           inputfile='/Users/fnespoli/run_dis/input_W7X_PMPI.txt'
                inputfile='/Users/fnespoli/run_dis/input_W7X.txt'
                inputfile='/Users/fnespoli/scripts/DIS/input.txt'
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

disp(['loading',folder2,fname0,'.mat'])
load([folder2,fname0,'.mat'])
   folder2=data{1}{18};
mi=1.67e-27;%kg
A=1;
rho_s=1.02*sqrt(T0*A(1))/(B0.*1e4); %m

if strcmp(material,'C')
    mI=12*1.67e-27;%kg
elseif strcmp(material,'B')
    mI=10.8*1.67e-27;%kg
elseif strcmp(material,'B4C')
    mI=55.255*1.67e-27;%kg
elseif strcmp(material,'BN')
    mI=24.82*1.67e-27;%kg
elseif strcmp(material,'Li')
    mI=6.94*1.67e-27;%kg
else
    disp([material, ' :material not recognized. Exiting'])
    if ~test exit; end
end
disp(['Computing neutral atom source, massrate=',num2str(massrate),'kg/s'])
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
