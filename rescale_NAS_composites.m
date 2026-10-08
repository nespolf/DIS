function [m1 m2 NAS1 NAS2]=rescale_NAS_composites(fnameBsource,material)
NAS=load(fnameBsource);
NAS1=NAS; NAS2=NAS;
 if strcmp(material,'BN')
     m1='B'; m2='N';
    % NAS1(:,2)=NAS(:,2).*10.8./24.82; %boron
    % NAS2(:,2)=NAS(:,2).*(24.82-10.8)./24.82;%nitrogen
     NAS1(:,2)=NAS(:,2); %boron
     NAS2(:,2)=NAS(:,2);%nitrogen
 elseif strcmp(material,'B4C')
     m1='B'; m2='C';
    % NAS1(:,2)=NAS(:,2).*4.*10.8./55.255;%boron
    % NAS2(:,2)=NAS(:,2).*(55.255-4.*10.8)./55.255;%carbon
      NAS1(:,2)=NAS(:,2).*4;%boron
      NAS2(:,2)=NAS(:,2);%carbon
    
     
 else disp('composite material not supported')
 end
 if exist('m1')
     fname1=[fnameBsource(1:end-4),'_',m1,'.txt'];
     fname2=[fnameBsource(1:end-4),'_',m2,'.txt'];
     dlmwrite(fname1,NAS1,'delimiter',' ')
     dlmwrite(fname2,NAS2,'delimiter',' ')
 end

end

