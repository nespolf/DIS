clear all
fname='~/dust_results_3/temp/JOREK_test_001.mat'
fname='~/dust_results_3/temp/160035_1.mat'

load(fname);

  rho_s=1.02*sqrt(2*T0)/(B0.*1e4);
  cs0=9.79e3*sqrt(T0);

  if material=='B'
      mI=10.8*1.67e-27;%kg
    rho_d=2080; % kg/m^3
  end
  
      m_d=4*pi/3*(a*rho_s).^3.*rho_d;
      Mu=Mu;
      for k=1:size(Mu,1)
      dMass(k,1)=0;
      dMass(k,2:size(Mu,2))=(Mu(k,2:end)-Mu(k,1:end-1)).*m_d;
      end
      
      for i=1:size(dMass,2)
          MassEj(i)=sum(dMass(1:i));
      end
      
      
      Massrate=-dMass./(dt.*rho_s./cs0);
      t=[1:length(Mu)].*(dt.*rho_s./cs0);
      AtomsPerSec=Massrate./mI;
      
      sourceEMC3=-dMass/mI./(t(end))*1.6e-19; %in [A];EMC3 is steady-state, so you have to assume that the source is spread evenly over the time of the trajectory?