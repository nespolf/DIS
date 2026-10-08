clear all
addpath /Users/fnespoli/scripts
addpath /Users/fnespoli/scripts/dust_dyn
cd /Users/fnespoli/TK3X
for i=1:10
    folder_list{i}= ['TCVlim_',num2str(20+i),'/Output/'];%
end

param = tk3x_load_param(folder_list{1});
mesh = tk3x_load_mesh(folder_list{1},param,0);
Bfield = tk3x_load_Bfield(folder_list{1},mesh,param.R0,0);
fields_list = {'N','Gammai'};%,'W'}%,'Ti','Te'}%,'Gammai','Gammae'}%,'Gammai'}%,'Te','Ti'}%,'PHI','Gammai','Gammae'};%,'W','Gammai','Te','Ti'};
Ntload_force = -1;plotbool =0; dtmin = 0;

            Nt0=0;
            dtload=4; 
            dtload=3000;
            [param,mesh,metric,Bfield,tt,fields] = tk3x_load_all(folder_list,fields_list,Nt0,dtload,plotbool,Ntload_force);

%%
for izone=1:length(mesh)
    clear N M deltaM deltaN
    N=squeeze(mean(mean(fields(izone).N,4),1));
    M=squeeze(mean(mean(fields(izone).Mi,4),1));
    NM=squeeze(mean(mean(fields(izone).Mi.*fields(izone).N,4),1));
   for it=1:size(fields(izone).N,1)
       for k=1:mesh(1).Nphi
       deltaN(it,:,:,k)=squeeze(fields(izone).N(it,:,:,k))-N;
       deltaM(it,:,:,k)=squeeze(fields(izone).Mi(it,:,:,k))-M;
       deltaNM(it,:,:,k)=squeeze(fields(izone).Mi(it,:,:,k).*fields(izone).N(it,:,:,k))-NM;
       end       
   end
   turb_flux(izone).phase=acos( mean(mean(deltaN.*deltaM,4),1)./sqrt(mean(mean(deltaN.^2,4),1))./sqrt(mean(mean(deltaM.^2,4),1))); 
end


 for i=1:length(mesh)
  phase(i).data=squeeze(turb_flux(i).phase)./pi;

 end
 h=figure;
 n =400;% get(gcf,'Number');
 output=tk3x_pol2D(phase,mesh,n,'\Gamma_\Psi'); %caxis([-1 1].*0.1)

mean(mean(squeeze(turb_flux(2).phase)./pi))