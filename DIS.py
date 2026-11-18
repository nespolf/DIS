#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#initialization
#globals().clear()

#Python translation of the open source release September 2026
#Library needed
import sys
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import scipy.special as scsp
import scipy.interpolate as scint
import scipy.integrate as scinteg

#Function definition
def cp_func(material, td, meltT, boilT):
   """
   Calculates the volumic mass, the heat capacity and the ablation enthalpy of the dust for each iteration
   """
   if material=='B': #https://webbook.nist.gov/
       mmol = 10.81
       alpha_v = 3*6.4e-6 
       rho_d = 2340 #T. Lundstrom Journal of Alloys and Compunds 267 1998
       if td>=meltT:
          rho_d = 2080
          CP = [31.75, 0, 0, 0, 0]
          alpha_v = 0
       # rho_d=(2.468-5.31e-5.*td).*1000; %F. Millot International Journal of Thermophysics, Vol. 23, No. 5, September 2002 ( 2002)
       elif td >= 1800:
          CP = [25.12664, 1.975493, 0.338395, -0.040032, -2.635578]
       else:
          CP = [10.18574, 29.24415, -18.02137, 4.212326, -0.550999]

   elif material=='W':
       mmol = 0.18384
       Rgp = 8.314
       alpha_v = 5E-6
       rho_d = 19300/(1+alpha_v*(td-298))**3
       cpT = lambda x: (21.452+1.24343e-2*x-1.23607e-5*x**2+7.86998e-9*x**3-2.27048e-12*x**4+2.67677e-16*x**5)/mmol
       cp = cpT(td)
       DeltacpT = lambda x: 3/2*Rgp/mmol-(21.452+1.24343e-2*x-1.23607e-5*x**2+7.86998e-9*x**3- 2.27048e-12*x**4+ 2.67677e-16*x**5)/mmol
       hsub, errhsub = scinteg.quad(DeltacpT, 298, td)
       hsub = hsub + 4670000
       if td>=meltT:
          alpha_v = 70E-6
          rho_d = 16300/(1+alpha_v*(td-meltT))**3
          cp = 291
          hsub = 4480000 + (3/2*Rgp/mmol-cp)*(boilT-td)
  
   elif material=='B4C': #https://webbook.nist.gov/
       mmol = 55.255; rho_d = 2520; #solid, cannot find liquid state
       alpha_v = 3*5.65e-6 #G.VTsagareishvili  Journal of the Less Common Metals 117 1986
       rho_d = rho_d/(1+alpha_v*(td-298))
       if td>=meltT:
          CP = [136, 0, 0, 0, 0]
       else:
          CP = [95.99853, 23.16513, -0.409604, 0.081414, -4.395208]
    
   elif material=='Li': #https://webbook.nist.gov/
       mmol = 6.94;
       alpha_v = 3*46.8e-6 #THERMAL EXPANSION OF LITHIUM, 77° TO 300° K.  WB Pearson · 1954
       alpha_v = 1.968e-4-6.724/td**2+1.413e5/td**4;# an evaluation of some thermodynamic and transport properties of solid and liquid lithium over the tempreature rnage 200-1700K R.K. Williams 1988
       if td>=meltT:
          rho_d = 512
          alpha_v = 1.01e-4/(0.5584-1.01e-4*td)
          CP = [29, 0, 0, 0, 0]
       else:
          rho_d = 534
          CP = [169.5520, -882.7110, 1977.438, -1487.312, -1.609635]
    
   elif material=='BN': #https://webbook.nist.gov/
       mmol = 24.82
       rho_d = 2100 #hexagonal, 3450 for c-BN
       rho_d = 2250 #from Momentive specs
       alpha_par = 38e-6 #hexagonal, M.E. Levinshtein, S.L. Rumyantsev, M.S. Shur Properties of Advanced Semiconductor Materials, 2001
       alpha_perp = -2.7e-6
       alpha_v = alpha_par+2*alpha_perp
    
       if td<1100:
          CP = [-1.626426, 101.9817, -81.52482, 25.65093, -0.223052]
       else:
          CP = [78.14750, -23.33367, 6.844396, -0.699753, -17.04034]
    
   elif material=='C': #Butland JNM 1973
       mmol = 12
       alpha_z = 27e-6+3.05e-9*(td-273) #thermal expansion coefficient
       alpha_b = -1.5e-6*(td<=523) + (1e-8*td-6.73e-6)*(td<=673)*(td>523)+ (0.45e-8*td-3.02e-6)*(td<=873)*(td>673) +(0.9e-6)*(td<=1073)*(td>873)+(2.5e-10*td+0.63e-6)*(td<=1273)*(td>1073)+(0.95e-6)*(td>1273)
       alpha_v = alpha_z+2*alpha_b
       rho_d = 2267 #graphite kg/m^3
       CP = [0.539, 9.11e-6, -90.27, -43449, 1.59e7, -1.4437e9] #Butland JNM 1973 cal/g/K
       CP = CP*4.1868*1000 #J/K/kg

   else:
       print('Temperature varying Cp for material '+material+' not implemented yet' )

   if material=='C':
       cpT = td
       cp = (CP[1] + CP[2]*cpT + CP[3]*(cpT)**-1 + CP[4]*(cpT)**-2 + CP[5]*(cpT)**-3 +CP[6]*(cpT)**-4) #function for boron
   elif material!='W':
       hsub = 508*1e3/10.811e-3
       cpT = td/1000
       cp = (CP[1] + (CP[2]*cpT) + (CP[3]*(cpT)**2) + (CP[4]*(cpT)**3) + (CP[5]/(cpT**2)))/(mmol*1e-3) #function for boron

   return rho_d, cp, hsub

def Coul_log(chi, tau, a, u, lD):
   """
   Calculates the Coulomb logarithm for drag forces
   """
   b0 = a*chi/tau/(3+2*u**2)
   ls2 = lD**2/(1+3/tau/(3+2*u**2))
   etafit = 1+a/np.sqrt(ls2)*(1+np.sqrt(1/(6*tau)))
   lnLambda = 0.5*np.log((b0**2+etafit**2*ls2)/(b0**2+a**2))
   return lnLambda

def Pvap(Td, Pvap_A, Pvap_B, material):
   """
   Calculates the vapor pressure of tungstene at dust temperature 
   """
   #Vapor pressure from Clausius-Clapeyron
   Pvap = 10**(Pvap_A/Td+Pvap_B) #[Pa] 

   #Vapor pressure from tabulation
   #if material=='W':
   #   M = np.genfromtxt('/home/AI289520/DIS_CODE/DIS/Pvap_W.csv')
    #  mindiffTd = (np.abs(M[:, 0]-Td)).argmin()
    #  Pvap = M[mindiffTd, 1]*1e5
   return Pvap

def vapor_shielding(Te, hsub, mI, kb, qe, ad, epsilon0, me):
   """
   Calculates the maximal ablation rate of the dust
   """
   #constants
   ZW = 15
   L0 = 1e-31
   RM = 2.5
   Eev = hsub*mI
   Ted = 1*qe
   Kew = 1e-13
   IW = 7.86*qe
   ne = 1E19

   #useful quantities
   Lambda = np.log(np.sqrt(2)*6*np.pi*(epsilon0*Te)**(3/2)/(ZW**2*np.sqrt(ne)*qe**2))
   kappae = 16*np.sqrt(2)*epsilon0**2*(Te)**(5/2)/(np.pi**(3/2)*ZW**2*qe**(3/2)*Lambda*np.sqrt(me))
   alphaL = kappae*L0*mI/(Eev**2)
   Gammak = 2*np.pi*ad**2*kappae*qe*Te/(Eev*np.sqrt(ad*RM))
   Gamma0 = 4*np.pi*ad*(Ted)**3/(me*Kew*IW*Eev)

   #Particule flux
   Gamma_vs = Gamma0**(3/5)*Gammak**(2/5)/alphaL**(1/5)*mI
   return Gamma_vs

def q_th(T, x, Wf, Te, me, kb, h, qe):
   """
   Calculates the heat flux by thermionic emission
   """
   #Te in eV
   Ithp = 16*np.pi*me*(kb*T)**2/h**3*np.exp(-Wf*qe/kb/T)

   if x>0:
      q_th = (2*T*kb/qe)*Ithp/4
   elif x<=0:
      q_th = (T*kb/qe)*Ithp/4*(2-2*x*Te*qe/kb/T+(x*Te*qe/kb/T)**2)*np.exp(x*Te*qe/kb/T)
   return q_th
   
def q_bb(T, Tw, epsilon_d, sigma_SB, qe):
   """
   Calculates the heat flux by grey body emission
   """
   return sigma_SB*epsilon_d*(T**4-Tw**4)/qe # last factor to convert to eV/m^2/s

def q_sub(T, Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te, epsilon0, me, m_d, cp, material):
   """
   Calculates the heat flux due to dust ablation
   """
   q_sub = Pvap(T, Pvap_A, Pvap_B, material)/np.sqrt(2*np.pi*mI*kb*T)*hsub/qe*(mI)

   if False: #q_sub>vapor_shielding(Te, hsub, mI, kb, qe, ad, epsilon0, me)*hsub/qe:
      q_sub = vapor_shielding(Te, hsub, mI, kb, qe, ad, epsilon0, me)*hsub/qe/m_d/cp
   return q_sub
   
def zeta_drag_n(u, v):
   """
   Calculates the neutral drag force
   """
   #from [Smirnov PPCF 2007]
   if v!=0:
      u = u/v
      zeta = 1/u*(1/np.sqrt(np.pi)*(u+1/2/u)*np.exp(-u**2)+(1+u**2-1/4/u**2)*scsp.erf(u))
   else:
      zeta = 0
   if np.isnan(zeta):
      zeta = 0
   return zeta
   
def zeta_drag(x, chi, tau, lnLambda):
   """
   Calculates the ion drag force 
   """
   # from [Smirnov PPCF 2007]
   xxp = x**2 + chi/tau
   xxm = x**2 - chi/tau
   xp = x + np.sqrt(np.abs(chi)/tau)
   xm = x - np.sqrt(np.abs(chi)/tau)
   if chi>=0:
      zeta_coll = 1/(2*x**2)*(1/np.sqrt(np.pi)*(1+2*xxp)*np.exp(-x**2)+x*(1+2*xxp-1/(2*x**2)*(1-2*xxm))*scsp.erf(x))
   else:
      zeta_coll = 1/(4*x**2)*(1/np.sqrt(np.pi)*((1+2*x**2+(1-2*x**2)/x*np.sqrt(-chi/tau))*np.exp(-xp**2) + (1+2*x**2-(1-2*x**2)/x*np.sqrt(-chi/tau))*np.exp(-xm**2)) + x*(1+2*xxp-1/(2*x**2)*(1-2*xxm))*(scsp.erf(xp)+scsp.erf(xm))) 
   zeta_orb = 2*lnLambda*(chi/tau)**2*(scsp.erf(x)-2*x/np.sqrt(np.pi)*np.exp(-x**2))/(2*x**3)
   zeta = zeta_coll + zeta_orb
   if np.isnan(zeta):
      zeta = 0
   return zeta
   
def eqTemp(y, t, ktemp, qred, x, sw_th, Te, Wf, Pvap_A, Pvap_B, hsub, mI, epsilon_d, Tw, kb, qe, me, h, sigma_SB, ad, epsilon0, m_d, cp, material):
   """
   Function to solve temperature equation
   """
   #all heat fluxes are devided by 4piad^2, units of eV/m^2/s
   dydt = ktemp*(qred - q_sub(y, Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te, epsilon0, me, m_d, cp, material) - sw_th*q_th(y, x, Wf, Te, me, kb, h, qe) - q_bb(y, Tw, epsilon_d, sigma_SB, qe))
   return dydt
   
def eqMass(y, t, td, k1, Pvap_A, Pvap_B, mI, kb, Te, hsub, qe, ad, epsilon0, me, mu, m_d, material):
   """
   Function to solve mass equation
   """
   gamma0 = Pvap(td, Pvap_A, Pvap_B, material)/np.sqrt(2*np.pi*mI*kb*td)*mI
   gamma_vs = vapor_shielding(Te, hsub, mI, kb, qe, ad, epsilon0, me)

   if True: #gamma0<gamma_vs:
      dydt = - Pvap(td, Pvap_A, Pvap_B, material)/np.sqrt(2*np.pi*mI*kb*td)*k1*y**(2/3)
   else:
      dydt = - gamma_vs/m_d*rho_s/cs0
   return dydt
   
def eqMot(y, t, Ni, M, k, mu, Ti, chi, lnLambda, Te, E, kE, A, Zi, Nn, Tn, Mn, mami, qw, fcent1, fcent2, fg):
   """
   Function to solve position equation
   """
   dydt = np.zeros(len(y))
   mid = int(len(y)/2)
   tot_l = int(len(y))
   dydt[:mid] = y[mid:tot_l]
   F = 0
   
   for z in range(len(Zi)):
      if Ni[z]>0:
         F = F + zeta_drag(np.abs(M-y[mid:tot_l])/np.sqrt(Ti[z]), Zi[z]*chi, Ti[z]/Te, lnLambda[z])*np.pi*mu**(2/3)*k*A[z]*Ni[z]*(M-y[mid:tot_l])*np.sqrt(Ti[z]*A[0]/A[z])
      else:
         pass
   if y[:mid]!=0:
      dydt[mid:tot_l] = 1/mu*F + 1/mu*zeta_drag_n(abs(Mn-y[mid:tot_l]), np.sqrt(Tn/mami))*np.pi*mu**(2/3)*k*mami*Nn*(Mn-y[mid:tot_l])*np.sqrt(Tn/mami) + 1/mu*(-chi*Te*mu**(1/3)*E*kE) + fcent1 + fcent2/y[:mid] + fg
   else:
      dydt[mid:tot_l] = 1/mu*F + 1/mu*zeta_drag_n(abs(Mn-y[mid:tot_l]), np.sqrt(Tn/mami))*np.pi*mu**(2/3)*k*mami*Nn*(Mn-y[mid:tot_l])*np.sqrt(Tn/mami) + 1/mu*(-chi*Te*mu**(1/3)*E*kE) + fcent1 + fg
   #qw/(4*pi*epsilon0)
   return dydt
   
def wall_reflection():
   """
   Calculates the new position of a dust when it hits a wall
   """
   return
   
def reassamble_dust_results(folder1, folder2, fname0, rd, A):
   """
   Reassemble all the dust trajectories in one file and compute the neutral atom source
   """
   print('Reassembling dust trajectories...')
   #Determine the shapes of the array that contain all the data
   fname1 = [fich for fich in os.listdir(folder1) if fich.startswith(fname0)]
   ktot = len(fname1)
   lll = np.zeros(ktot, dtype = int)
   for k in range(ktot): 
      datak = scio.loadmat(folder1+fname1[k])
      lll[k] = np.shape(datak['R'])[1]
      if k==0:
         Nd = datak['Nd']
   
   #Creation of the arrays
   R1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan 
   Z1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   phi1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   vphi1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   vZ1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   vR1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   Mu1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   Nd1 = np.reshape(np.zeros(ktot*max(lll)*np.shape(Nd)[1]), (ktot,max(lll),np.shape(Nd)[1]))*np.nan 
   Td1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   Ted1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   Tau1 = np.reshape(np.zeros(ktot*max(lll)*np.shape(Nd)[1]), (ktot,max(lll),np.shape(Nd)[1]))*np.nan
   Chi1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   Md1 = np.reshape(np.zeros(ktot*max(lll)*3), (ktot,max(lll),3))*np.nan
   Ed1 = np.reshape(np.zeros(ktot*max(lll)*3), (ktot,max(lll),3))*np.nan
   t1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   dead1 = np.array(range(1, ktot+1))*np.nan
   Rho_d1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   Cp1 = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   dMu = np.reshape(np.zeros(ktot*max(lll)), (ktot, max(lll)))*np.nan
   
   #Giving the values for each dust trajectory
   for k in range(ktot):
      datak = scio.loadmat(folder1+fname1[k])
      R1[k, :lll[k]] = datak['R']; Z1[k, :lll[k]] = datak['Z']; phi1[k, :lll[k]] = datak['phi'];
      vR1[k, :lll[k]] = datak['vR']; vZ1[k, :lll[k]] = datak['vZ']; vphi1[k, :lll[k]] = datak['vphi'];
      Mu1[k, :lll[k]] = datak['Mu']; Td1[k, :lll[k]] = datak['Td']; Chi1[k, :lll[k]] = datak['Chi'];
      Nd1[k, :lll[k], :] = datak['Nd']; Ted1[k, :lll[k]] = datak['Ted']; Tau1[k, :lll[k], :] = datak['Tau'];
      Md1[k, :lll[k], :] = datak['Md']; Ed1[k, :lll[k], :] = datak['Ed']; t1[k, :lll[k]] = datak['t'];
      dead1[k] = datak['dead']; Rho_d1[k, :lll[k]] = datak['Rho_d']; Cp1[k, :lll[k]] = datak['Cp'];

   R = R1; Z = Z1; phi = phi1;
   vR = vR1; vZ = vZ1; vphi = vphi1;
   Mu = Mu1; Td = Td1; Chi = Chi1;
   Nd = Nd1; Md = Md1; Ed = Ed1;
   Ted = Ted1; Tau = Tau1; t = t1;
   dead = dead1; Cp = Cp1; Rho_d = Rho_d1;
   rho_s = datak['rho_s']
   cs0 = datak['cs0']
   
   if np.sum(Ed)==0:
      Ed = 0
      print('No electric field, zero everywhere.')

   print('Mass lost at the outside of the plasma: '+f"{(1-np.abs(np.nansum(np.diff(Mu))))*100:.4f}"+'%') #_end
   del R1, Z1, phi1, vphi1, vZ1, vR1, Mu1, Nd1, Md1, Td1, Ted1, Tau1, Chi1, Ed1, t1, dead1, Cp1, Rho_d1, k, ktot, datak#, iii, aaa
   
   #Computing the neutral atom source for all each dust
   NA = 6.02e23
   dt = np.zeros(np.shape(Mu)[1])
   if np.shape(Mu)[1]==1 and np.shape(Mu)[0]>1:
      Mu = np.transpose(Mu)
   for i in range(np.shape(Mu)[0]):
      for it in range(np.shape(Mu)[1]):
         dMu[i, it] = Mu[i, it] - Mu[i, it-1]
         dt[it] = (t[i, it] - t[i, it-1])*rho_s/cs0
   dMu = 1
   NAS = - 4/3*np.pi*rd**3*Rho_d*dMu/A*1000*NA/dt
   
   fnameout = folder2+fname0+'.mat'
   scio.savemat(fnameout, {'NKU' : NKU, 'n0' : n0, 'T0' : T0, 'B0' : B0, 'cs0' : cs0, 'rho_s' : rho_s, 'a' : a, 'x0' : x0, 'v0' : v0, 'R' : R, 'vR' : vR, 'Z' : Z, 'vZ' : vZ, 'phi' : phi, 'vphi' : vphi, 'Chi' : Chi, 'Qw' : Qw, 'Rd' : Rd,'Tau' : Tau, 'Td' : Td, 'Mu' : Mu, 'ITHd' : ITHd, 'IRd' : IRd, 'IZd' : IZd, 'IPHId' : IPHId, 'Nd' : Nd, 'Ted' : Ted, 'Md' : Md, 'Ed' : Ed, 'material' : material, 't' : t, 'dt' : dt, 'test' : test, 'dead' : dead, 'Psimelt' : Psimelt, 'Psiboil' : Psiboil, 'Nnd' : Nnd, 'Tnd' : Tnd, 'Rho_d' : Rho_d, 'Cp' : Cp, 'NAS' : NAS})
   print('File saved in '+fnameout)
   return NAS
   
def outputdata_HFPS():
   """
   Generates a file that computes the relevant output for HFPS
   """
   return 

#Main code
if __name__ == "__main__":
	print('*******************************************')
	print('*******************************************')
	print('DIS starts running')
	print('-------------------------------------------')
	print('-------------------------------------------')
	
	#Version de python exécutée
	print('The executed python version is :')
	print(sys.executable)
	print('-------------------------------------------')

	test = False
        
	if test:
           mainfolder = '/home/AI289520/DIS_CODE/DIS/' 
           inputfile = 'input_DIS.txt'
	else: 
           mainfolder = '/home/AI289520/DIS_CODE/DIS/' 
           inputfile ='input.txt'
           print('Uploading the input data...')
           print('The input file is: '+mainfolder+inputfile)
           
        #Give access to all the functions and files in the mainfolder
	sys.path.append(mainfolder)
	
	#Launching the clock
	times = []
	t0 = time.time()
	times.append(t0)
        
        #Getting all the data from the input file
	with open(mainfolder+inputfile, 'r', encoding='utf-8') as f:
	   data = f.readlines()
	   
	rd = float(data[1].strip())
	material = data[3].strip()
	x0i = [float(x) for x in data[5].strip().split()]
	v0i = [float(x) for x in data[7].strip().split()]
	NKU = float(data[9].strip())
	phi_or_t = [float(x) for x in data[11].strip().split()]
	if len(phi_or_t)>1:
	   alpha_0 = phi_or_t[2] 
	   phi_or_t = phi_or_t[1]
	   print('Randomized velocity with std of '+str(alpha_0)+' deg')
	else:
	   phi_or_t = phi_or_t[0]
	   if phi_or_t==4:
	      alpha_0 = 10
	fnamepbg = data[13].strip()
	folder1 = data[15].strip()
	folder2 = data[17].strip()
	fname0 = data[19].strip()
	xi = round(np.sqrt(x0i[0]**2+x0i[1]**2+x0i[2]**2)*1E2)
	vi = round(np.sqrt(v0i[0]**2+v0i[1]**2+v0i[2]**2)*1E2)
	fname0 = fname0+'_'+material+'_'+str(rd)+'m_'+str(xi)+'m_'+str(vi)+'ms-1'
	fnamewall = data[21].strip()
	if data[23].strip()=='NO':
	   massrate = []
	else:
	   massrate = float(data[23].strip())
	sw = [float(x) for x in data[25].strip().split()]
	sw_neu = sw[0] 
	sw_imp = sw[1]
	mainion = data[27].strip()

	#Upload plasma background
	print('Uploading plasma background...')
	print('The background file is: '+fnamepbg)
	print('-------------------------------------------')
	pbg = scio.loadmat(fnamepbg)

	B0 = pbg['B0'][0][0]
	T0 = pbg['T0'][0][0]
	cs0 = pbg['cs0'][0][0]
	n0 = pbg['n0'][0][0]
	fieldER = pbg['fieldER']
	fieldEZ = pbg['fieldEZ']
	fieldEphi = pbg['fieldEphi']
	fieldNe = pbg['fieldN']/n0
	fieldNi = pbg['fieldN']/n0
	fieldTe = pbg['fieldTe']/T0
	fieldTi = pbg['fieldTi']/T0
	fieldmR = pbg['fieldmR']
	fieldmZ = pbg['fieldmZ']
	fieldmphi = pbg['fieldmphi']
	param = pbg['param']
	mesh = pbg['mesh'][0]
	meshR, meshZ, meshphi = mesh[0]

	if np.shape(fieldNi)[0]==1:
	   print('No impurities provided in the background plasma. Impurities excluded.')
	elif sw_imp!=1:
	   fieldNi = fieldNi[0,:,:,:] #to not load impurities
	   print('Impurities excluded.')
	   
	#Getting the neutral background
	try:
	   fieldNn=fieldNn/n0
	   fieldTn=fieldTn/T0
	   if sw_neu!=1:
	      del fieldNn, fieldTn
	      print('Neutrals excluded.')
	except NameError:
	   print('No neutral density  provided in background plasma. Neutrals excluded.')
	
	#Get the nuclei parameters for ion
	mami = 1 #atom/ion mass ratio
	if np.shape(fieldNi)[0]==1: #H or D
	   A = [1]
	   Zi = [1]
	elif np.shape(fieldNi)[0]==2: #He
	   A = [4, 4]
	   Zi = [1, 2]
	elif np.shape(fieldNi)[0]==7: #H or D+C
	   A = [1, 12, 12, 12, 12, 12, 12]
	   Zi = [1, 1, 2, 3, 4, 5, 6]
	if mainion=='H':
	   A[0] = 1
	elif mainion=='D':
	   A[0] = 2
	else:
	   print('Main ion not recongnized/implemented')

	#Access the walls coordinates
	if fnamewall=='NO':
	   print('No wall file specified, reflections with wall ignored.')
	try:
	   wall = scio.loadmat(fnamewall) #wall units in cm
	   Rwall = wall['Rwall'][:, 0]/100
	   Zwall = wall['Rwall'][:, 1]/100
	except FileNotFoundError:
	   pass

	#intialization: setting the values of the physical constants
	epsilon0 = 8.8542e-12
	mi = A[0]*1.67e-27
	me = 9.1094e-31
	Tw = 300 #340K pour le divertor activement refroidit de WEST ???
	qe = 1.6022e-19
	B0 = 1
	g_T = 9.81
	sigma_SB = 5.67e-8
	h = 6.6261e-34
	kb = 1.3807e-23
	mime = np.sqrt(mi/me)
	cs0 = 9.79e3*np.sqrt(2*T0/A[0])
	rho_s = 1.02*np.sqrt(T0*A[0])/(B0*1e4)
	E_rec = 13.6 #eV, ionization potentia energy (H) for recombination on dust surface + H2 molecule formation
	E_mol = 1.1 #eV, energy of molecular recombination
 
	#Adimensionalisation of the mesh
	meshR = meshR/rho_s
	meshZ = meshZ/rho_s

	if len(meshphi)>1:
	   print('Case not taken into account for the moment')
	   sys.exit()
	   #mesh(1).phinodes=mesh(1).phicells-mean(diff(mesh(1).phicells))/2;
	   #mesh(1).phinodes(end+1)=mesh(1).phinodes(end)+diff(mesh(1).phinodes(1:2));
	else:
	   meshphinodes = [0, 1]

	if max(meshphi)>2*np.pi:
	   print('Case not taken into account for the moment')
	   #meshphinodes = mesh(1).phinodes*pi/180;
	   #meshphinodes = mesh(1).phicells*pi/180;

	#Speed-position initial conditions
	x0i = [x0i[0]/rho_s, x0i[1]/rho_s, x0i[2]]
	X0 = x0i[0]*np.cos(x0i[2])
	Y0 = x0i[0]*np.sin(x0i[2])
	Z0 = x0i[1]
	v0i = v0i/cs0
	v00 = np.sqrt(v0i[0]**2+v0i[1]**2+v0i[2]**2)
 
	#Step time normalisation
	dt = round(5e-6/(rho_s/cs0)) # dt=5e-6 s, normalized
	t_f = np.linspace(0, dt, np.shape(fieldNe)[1]+1)

	#in this code M is Vp/cs0=Vp/cs*(cs/cs0)
	fieldmR = fieldmR*np.sqrt(fieldTe)
	fieldmZ = fieldmZ*np.sqrt(fieldTe)
	fieldmphi = fieldmphi*np.sqrt(fieldTe)
	
	try:
	   fieldER = fieldER/T0*rho_s
	   fieldEZ = fieldEZ/T0*rho_s
	   fieldEphi = fieldEphi/T0*rho_s
	except NameError:
	   pass

	#initialization: giving the right thermodynamical values for each material
	if material=='C':
	   Pvap_A = -40181 
	   Pvap_B = 14.8 
	   mI = 12*1.67e-27
	   rho_d = 2267
	   cp = 1709
	   Wf = 4.6
	   hsub = 2.97e7
	   epsilon_d = 0.75
	   b_see = 0.4
	   boilT = 3925
	   meltT = 1e6
	elif material=='B':
	   Pvap_A=-28238
	   mI=10.8*1.67e-27
	   rho_d = 2080
	   rho_d = 2340
	   cp = 11.087/10.811e-3
	   Wf = 4.45
	   hsub = 508*1e3/10.811e-3
	   hmel = 50.2*1e3/10.811e-3
	   epsilon_d = 0.27
	   epsilon_d = 0.8
	   b_see = 0.5
	   boilT = 4200
	   meltT = 2349
	elif material=='W':
	   Pvap_A = -45100
	   Pvap_B = 12.85
	   mI = 183.84*1.67e-27
	   rho_d = 19300 
	   cp = 132
	   Wf = 4.55
	   heb = 4480000
	   hsub = 4670000
	   hmel = 192000
	   epsilon_d = 0.35
	   b_see = 0.8
	   boilT = 5828
	   meltT = 3695
	elif material=='B4C':
	   Pvap_A = -31690 
	   Pvap_B = 13.42
	   mI = 55.255*1.67e-27
	   rho_d = 2520
	   cp = 100/55.255e-3
	   Wf = (5.15+5.46)/2
	   hsub = 138.*4.184*1e3/55.255e-3
	   hmel = (1350+2030)/2.*1e3
	   epsilon_d = 0.93
	   b_see = 0.5
	   boilT = 3770
	   meltT = 3063
	elif material=='BN':
	   Pvap_A = -25710
	   mI = 24.82*1.67e-27
	   rho_d = 2250
	   cp = 45/24.82e-3
	   Wf = 4.5
	   hmel = 81.*1e3/24.82e-3
	   hsub = hmel
	   epsilon_d = 0.85
	   b_see = 0.5
	   boilT = 3246
	   meltT = 1e6
	elif material=='Li':
	   Pvap_A = -8042 
	   Pvap_B = 10.02
	   mI = 6.94*1.67e-27
	   rho_d = 534
	   cp = 24.860/6.94e-3
	   cp = 27/6.94e-3
	   Wf = 2.93
	   hsub = 136*1e3/6.94e-3
	   hmel = 3.0*1e3/6.94e-3
	   epsilon_d = 0.04
	   b_see = 1.3
	   boilT = 1603
	   meltT = 453.65
	else:
	   print(material+' :material not recognized. Exiting.')
	   sys.exit()

	#Coefficients nedeed for SEE
	C_see = np.genfromtxt(mainfolder+'YoungDekkerFit_'+material+'.txt')
    
	#Normalization/switches
	a = rd/rho_s
	m_d = 4*np.pi/3*(a*rho_s)**3*rho_d
	kg = rho_s/cs0**2*g_T
	sw_th = 1
	sw_see = 1
	sw_BS = 0
	sw_REC = 0
    
	#Set of adimensionned dust potential
	x = np.linspace(-8, 3, 5000)
	HeavisideX = np.heaviside(x, 1)
	HeavisidemX = np.heaviside(-x, 1)
    
	#Averaging 3d turbulence
	fields_list = ['fieldNe','fieldNi','fieldTe','fieldTi','fieldER','fieldEZ','fieldEphi','fieldmR','fieldmZ','fieldmphi']
	if NKU==np.inf:
	   print('averaging 3D fields to 2D...')
	   for i in range(len(fields_list)):
	      try:
	         eval(fields_list[i]+'=(np.nanmean(np.nanmean('+fields_list[i]+',axis = 1), axis = 4))')
	      except NameError:
	         pass
	elif abs(NKU)>0:
	   print('averaging 3D fields to alter Ku number...')
	   for i in range(len(fields_list)):
	      try:
	         eval(f1[i]+'=(np.nanmean(np.nanmean('+fields_list[i]+',axis = 1), axis = 4))')
	         eval(fields_list[i]+'=('+fields_list[i]+'-f1)./(NKU+1)+f1')
	         del f1
	      except NameError:
	         pass
	else:
	   print('3D turbulence averaging not implemented yet.')

	#Ending initialization
	t1 = time.time()
	times.append(t1)
	tinit = t1 - t0
	print('-------------------------------------------')
	print('Time for initialization '+f"{tinit:.2f}"+' s')
	print('-------------------------------------------')
	print('-------------------------------------------')

	#La danse commence...
	print('Computing dust dynamics...')
  
	#initialisation of the number of trajectories with phi_or_t : 
	if phi_or_t==1: #change initial toroidal position for 2 tests and 4 for testing several directions
	   kkktot = len(meshphi)
	elif phi_or_t==2 and np.shape(fieldNe)[1]==1:
	   kkktot = 1
	elif phi_or_t==2 or phi_or_t==3:  #change initial time, or both, randomly
	   kkktot = 1000
	elif phi_or_t==4:  #random initial velocity direction
	   kkktot = 500
	elif phi_or_t==5: #scan R,phi position where plasma is actually there        
	   n1 = np.squeeze(np.nansum(np.squeeze(fieldNe))[1])
	   N = 10 #resampling in R because too many points
	   ilist = [] 
	   jlist = []
	   delta = meshR[::N, 0] - meshR[i, 0]
	   for i in range(np.shape(n1)[0]):
	      for j in range(np.shape(n1)[1]):
	         if n1[i, j]>0 and any(k == meshR[i, 0] for k in meshR[::N, 0]):
	            ilist.append(i)
	            jlist.append(j)
	   kkktot = len(ilist)
	   del N, n1

	if NKU==np.inf or test:
	   kkktot=1
    
	#solve for each dust
	for kkk in range(kkktot):
	   if NKU==np.inf or test:
	      kkk = 1
	   else:
	      fname1 = [fich for fich in os.listdir(folder1) if fich.startswith(fname0)]
	      if len(fname1)==0: 
	         kstart1 = 1
	      else:
	         nn = []
	         for l in range(len(fname1)):
	            num = fname1[l]
	            nn.append(int(num[-7:-4]))
	         klist = np.array(range(1, kkktot+1))
	         kleft = np.setdiff1d(klist, nn)
	         if len(kleft)==0:
	            print('The data files have already been computed in '+folder1)
	            print('-------------------------------------------')
	            break
	         else:
	            print('First missing file detected: '+str(kleft[0]))
	            kkk = kleft[0]
        
	   #Displaying the number of dust trajectory
	   fnameout = folder1 + fname0 + '_' + f"{kkk+1:03d}" + '.mat'
	   #scio.savemat(fnameout, {"kkk" : kkk+1}) place holder file
	   t2 = time.time()
	   times.append(t2)
	   print('Trajectory '+str(kkk+1)+'/'+str(kkktot))
     
	   #initialisation of the spatio-speed initial conditions with phi_or_t 
	   Nfact = 2*np.pi/(meshphinodes[-1]-meshphinodes[0])
	   if phi_or_t==1: #change the toroidal position
	      itf0 = 0
	      x0 = [x0i[0], x0i[1], x0i[2]+kkk*np.mean(np.diff(meshphi))]
	   elif phi_or_t==2: #change initial time
	      if len(t_f)==1:
	         tf0 = 1
	      else:
	         tf0 = 0 #np.random.rand(1)*(t_f[-1]-t_f[0])+t_f[0]
	         itf0 = (np.abs(t_f-tf0)).argmin()
	   elif phi_or_t==3: #both
	      tf0 = np.random.rand(1)*(t_f[-1]-t_f[0])+t_f[0]
	      itf0 = (np.abs(t_f-tf0)).argmin()
	      x0 = [x0i[0], x0i[1], meshphi[round(np.random.rand(1)*(len(meshphi)-1)+1)]]
	      print('itf0='+str(itf0)+', iphi0='+str(np.unravel_index((np.abs(meshphi-np.mod(x0[2], 2*np.pi/Nfact))).argmin(), np.shape(meshphi))))
	   elif phi_or_t==4: #randomize initial velocity direction
	      itf0 = 0
	      alpha = 180+alpha_0*np.random.randn(1) #normally distributed
	      beta = np.random.rand(1)*2*180-180 #uniformily distributed
	      v0 = np.zeros(3)
	      v0[0] = v00*np.sin(alpha*np.pi/180)*np.cos(beta*np.pi/180)
	      v0[1] = v00*cos(alpha*np.pi/180)
	      v0[2] = v00*sin(alpha*np.pi/180)*sin(beta*np.pi/180)     
	      X00 = X0 + 5e-3/rho_s*np.random.randn(1)
	      Y00 = Y0 + 5e-3/rho_s*np.random.randn(1)
	      Z00 = Z0
	      if v0[0]!=0:
	         X00 = X0
	         Y00 = Y0 + 5e-3/rho_s*np.random.randn(1)
	         Z00 = Z0 + 5e-3/rho_s*np.random.randn(1)
	         alpha = 180 + alpha_0*np.random.randn(1)
	         v0[0] = v00*cosd(alpha*np.pi/180)
	         v0[1] = v00*sind(alpha*np.pi/180)*cosd(beta*np.pi/180)
	         v0[2] = v00*sind(alpha*np.pi/180)*sind(beta*np.pi/180)
	         if v0i[1]!=0 or v0i[2]!=0:
	            aa = np.array([-1, 0, 0])
	            bb = v0i/np.sqrt(v0i[0]**2+v0i[1]**2+v0i[2]**2)
	            R = 2*np.outer(aa+bb, aa+bb)/np.outer(aa+bb, aa+bb)-np.identity(3)
	            v0 = R*v0
	            X00 = X0
	            Y00 = Y0
	            Z00 = Z0
	      x0 = [np.sqrt(X00**2+Y00**2), Z00, np.arctan(Y00,X00)]
	      print('x0 = ['+num2str(x0)+'], v0 = ['+num2str(v0)+']')
	   elif phi_or_t==5: #%scan radial and toroidal
	      itf0 = 0
	      x0 = [meshR(ilist[kkk]), x0i[1], meshphi(jlist[kkk])];
	      print('i_R = '+str(ilist[kkk])+'/'+str(len(meshR)),+' i_phi = '+str(jlist[kkk])+'/'+str(len(meshphi)))
        
	   #Clearing and re-initialisation of the variables were the dust dynamics is stored
	   try:
	      del t, R, vR, Z, vZ, phi, vphi, Chi, Qw, Rd, Tau, Td, Mu, ITHd, IRd, IZd, IPHId, Nd, Ted, Md, ERd, EZd, Ephid, Fi, FE, dead, MTOR, LNLAMBDA
	   except NameError:
	      pass
	   
	   #Initialisation of the vectors containing the computed data
	   Ntmax = 400000
	   R = np.zeros(Ntmax)*np.nan
	   t = np.zeros(Ntmax)*np.nan
	   vR = np.zeros(Ntmax)*np.nan
	   Z = np.zeros(Ntmax)*np.nan
	   vZ = np.zeros(Ntmax)*np.nan
	   phi = np.zeros(Ntmax)*np.nan
	   vphi = np.zeros(Ntmax)*np.nan
	   Chi = np.zeros(Ntmax)*np.nan
	   Qw = np.zeros(Ntmax)*np.nan
	   Rd = np.zeros(Ntmax)*np.nan
	   Tau = np.reshape(np.zeros(Ntmax*np.shape(fieldNi)[0]), (Ntmax, np.shape(fieldNi)[0]))*np.nan
	   Td = np.zeros(Ntmax)*np.nan
	   Mu = np.zeros(Ntmax)*np.nan
	   IPHId = np.zeros(Ntmax, dtype=int)
	   IRd = np.zeros(Ntmax, dtype=int)
	   ITHd = np.zeros(Ntmax, dtype=int)
	   IZd = np.zeros(Ntmax, dtype=int)
	   Nd = np.reshape(np.zeros(Ntmax*(np.shape(fieldNi)[0]+1)), (Ntmax, np.shape(fieldNi)[0]+1))*np.nan
	   Ted = np.zeros(Ntmax)*np.nan
	   Md = np.reshape(np.zeros(3*Ntmax), (Ntmax, 3))*np.nan
	   Ed = np.reshape(np.zeros(3*Ntmax), (Ntmax, 3))*np.nan
	   Psimelt = np.zeros(Ntmax)*np.nan
	   Psiboil = np.zeros(Ntmax)*np.nan
	   Nnd = np.zeros(Ntmax)*np.nan 
	   Tnd = np.zeros(Ntmax)*np.nan     
	   Rho_d = np.zeros(Ntmax)*np.nan  
	   Cp = np.zeros(Ntmax)*np.nan

	   #Initialization at t=0
	   solvedyn = True
	   it = 0
	   hasmelt = False
        
	   #Beginning of the dynamics resolution
	   while solvedyn:
	      if it==0:
	         print('Loop initialization...')
	         try:
	            x0
	         except NameError:
	            x0 = x0i
	         try:
	            v0
	         except NameError:
	            v0 = v0i
	         itf = itf0
	         R[it] = x0[0]
	         vR[it] = v0[0]
	         Z[it] = x0[1]
	         vZ[it] = v0[1]
	         phi[it] = x0[2]
	         vphi[it] = v0[2]
	         Chi[it] = 3 #-log(4*pi*me/mi*(1+5/3*Ti/Te))/2
	         Cp[it] = cp
	         Rho_d[it] = rho_d
	         Qw[it] = -1E-11 #- 4*np.pi*epsilon0*ad*(1+ad/lDtot)*Chi[it]*Te
	         Rd[it] = rd
	         Tau[it, :] = np.ones(len(A))
	         Td[it] = 300 
	         mu = 1
	         Mu[it] = mu
	         IPHId[it] = (np.abs(meshphi - x0[2])).argmin()
	         IRd[it] = (np.abs(meshR[:, 0] - x0[0])).argmin()
	         ITHd[it] = (np.abs(meshZ[0, :] - x0[1])).argmin()
	         IZd[it] = 0
	         Nd[it, :] = np.nan*np.ones(np.shape(fieldNi)[0]+1)
	         Ted[it] = fieldTe[0, itf, IRd[it], ITHd[it]]
	         Nnd[it] = np.nan
	         Tnd[it] = np.nan
	         Md[it] = [0, 0, 0]
	         Ed[it] = [0, 0, 0]
	         #Fi[it] = np.nan
	         #FE[it] = np.nan
	         Psimelt[it] = 0
	         psimelt = 0
	         Psiboil[it] = 0
	         psiboil = 0
	         t[it] = t_f[itf]
	         print('Resolving physics...')
            
	      #Taking the background values at particle position
	      ir = IRd[it]; ith = ITHd[it]; iz = IZd[it];
	      if min(meshphi)+max(meshphi)==0:
	         iphi = (np.abs(np.mod(phi[it]-np.pi/Nfact, 2*np.pi/Nfact) - np.pi/Nfact - meshphi)).argmin()
	      elif min(meshphi)>=0:
	         iphi = (np.abs(np.mod(phi[it], 2*np.pi/Nfact) - np.pi/Nfact - meshphi)).argmin()

	      if NKU!=np.inf:
	         iphi1 = iphi
	         itf = (np.abs(np.mod(t[it], t_f[-1]) - t_f)).argmin()
	      else:
	         itf = 0; iphi1 = 0;
            
	      #ITF[it] = itf;
	      Ni = fieldNi[:, itf, ir, ith]
	      Ne = fieldNe[0, itf, ir, ith]
	      if Ne<0:
	         dead = -1 
	         print('Negative plasma density')
	         break
	      Te = fieldTe[0, itf, ir, ith]
	      Ti = fieldTi[:, itf, ir, ith]
	      try:
	         Nn = fieldNn[0, itf, ir, ith]
	         Tn = fieldTn[0, itf, ir, ith]
	         Mn = 0
	      except NameError:
	         Nn = 0; Tn=0; Mn=0;

	      try:
	         ER = fieldER[0, itf, ir, ith]
	         EZ = fieldEZ[0, itf, ir, ith]
	         Ephi = fieldEphi[0, itf, ir, ith]
	      except NameError:
	         ER = 0; EZ = 0; Ephi = 0;
            
	      mR = fieldmR[iz, itf, ir, ith]
	      mZ = fieldmZ[iz, itf, ir, ith]
	      mphi = fieldmphi[iz, itf, ir, ith]
	      
	      #temperature dependent density and thermal capacity
	      rho_d, cp, hsub = cp_func(material, Td[it], meltT, boilT)
	      if np.isnan(Ne) or Te*T0<0.1:
	         mR = 0; mZ = 0; mphi = 0;  Nn = 0; Tn = 0; Ne = 1e13/n0; Te = 0.1/T0; Ti = Te; 
	         Ni = Ne*np.reshape(np.ones(len(Zi)), (1, len(Zi)))

	      if len(A)==1:
	         Ti = Ti*np.reshape(np.ones(len(Ni)), np.shape(Ni))
	      tau = Ti/Te
	      u = np.sqrt((mR-vR[it])**2+(mZ-vZ[it])**2+(mphi-vphi[it])**2)/np.sqrt(Ti)*np.sqrt(np.array(A)*A[0]**-1)
	      lD = np.sqrt(epsilon0*Te*T0/(qe*Ne*n0))
	      lDi = np.sqrt(epsilon0*Ti*T0/(qe*Ne*n0))
	      lDtot = np.sqrt(1/lD**2+1/lDi**2)
	      ad = (3/4/np.pi*m_d*mu/rho_d)**(1/3);
	      qw = -4*np.pi*epsilon0*ad*(1+ad/lDtot)*Chi[it]*Te
	      lnLambda = Coul_log(Chi[it], tau, ad, u, lD)
	      
	      #Normalization
	      ktemp = (4*np.pi/mu/m_d)**(1/3)*(3/rho_d)**(2/3)/cp*(rho_s/cs0)*qe
	      k = (3/(4*np.pi*rho_d))**(2/3)*mi*m_d**(-1/3)*n0*rho_s
	      k1 = (3/rho_d)**(2/3)*(4*np.pi/m_d)**(1/3)*mI*rho_s/cs0
	      kE = (3/4/np.pi/rho_d/m_d**2)**(1/3)*(T0/cs0)**2*(4*np.pi*epsilon0)
	      kL = (3/rho_d/m_d**2*(4*np.pi)**2)**(1/3)*(T0*B0*rho_s/cs0)*(epsilon0)
	      
	      #dust charge : all currents are divided by pi*a_d^2*e*cs0*n_e
	      T = Td[it]
	      
	      #electrons
	      rho_e = np.sqrt(me*qe*Te*T0)/(qe*B0)
	      rc = 0.3E-6
	      ksi = 1 #np.tanh((rho_e-ad)/rc)
	      Ie0 = ksi*np.sqrt(Te*4/np.pi)*mime+(1-ksi)*8*Te*T0/(B0*qe*np.pi*ad*cs0)*(1+Ti/Te)
	      Iep = - Ie0*np.exp(-x)
	      Ien = - Ie0*(1-x)
	      
	      #ions
	      Iip = 0; Iin = 0;
	      for z in range(len(Zi)):
	         Iip = Iip + Ni[z]/Ne*Zi[z]*np.sqrt(Ti[z]*A[0]/A[z])*(scsp.erf(u[z])*(u[z]+1/(2*u[z])+x*Zi[z]/(tau[z]*u[z]))+1/np.sqrt(np.pi)*np.exp(-u[z]**2))
	         um = u[z] - np.sqrt(Zi[z]*np.abs(x)/tau[z])
	         up = u[z] + np.sqrt(Zi[z]*np.abs(x)/tau[z])
	         Iin = Iin + Ni[z]/Ne*Zi[z]*np.sqrt(Ti[z]*A[0]/A[z])/4/u[z]*((scsp.erf(up[z]) + scsp.erf(um[z]))*(1+2*u[z]**2+2*x*Zi[z]/tau[z])+2/np.sqrt(np.pi)*(np.exp(-(um[z])**2)*(up[z])+np.exp(-(up[z])**2)*(um[z])))

	      #thermionic emission
	      Ithp = sw_th*16*np.pi*me*(kb*T)**2/h**3/cs0/Ne/n0*np.exp(-Wf*qe/kb/T)
	      Ithn = np.zeros(5000)
	      Ithn[:int(8/11*5000)] = Ithp*(1-x[:int(8/11*5000)]*(Te*T0*qe)/(T*kb))*np.exp(x[:int(8/11*5000)]*(Te*T0*qe)/(T*kb))
         
	      #SEE emission
	      F_f = 10**(C_see[0] + C_see[1]*np.log10(Te*T0) + C_see[2]*np.log10(Te*T0)**2 + C_see[3]*np.log10(Te*T0)**3 + C_see[4]*np.log10(Te*T0)**4)
	      dSEEp = sw_see*2/(2-b_see)*F_f
	      dSEEn = sw_see*2/(2-b_see)*F_f*(1-3*x*Te*T0/Wf)/(1-x)/(1-x*Te*T0/Wf)**3*np.exp(-x) #Autricque These
         
	      #Electron backscattering
	      fch = 0.60*(1-1.86*(Te*T0)**(-0.23)*np.exp(-0.0033*Te*T0))
	      IEBS = - sw_BS*Ie0*np.exp(-x)*fch
   
	      #Total current balance
	      eqtot = (Iep*(1-dSEEp)-IEBS+Iip+Ithp)*HeavisideX + (Ien*(1-dSEEn)-IEBS+Iin+Ithn)*HeavisidemX
 
 
	      #dust temperature : all fluxes are divided by 4pia^2
	      qip = 0
	      qin = 0
	      
	      #sum contribution for various ions
	      for z in range(len(A)):
	         qip = qip + Ni[z]*np.sqrt(Ti[z]*A[0]/A[z])*Ti[z]*1/4*(2/np.sqrt(np.pi)*(5+2*u[z]**2+2*x*Zi[z]/tau[z])*np.exp(-u[z]**2)+(3/u[z]+12*u[z]+4*u[z]**3+2*x*Zi[z]/tau[z]*(1/u[z]+2*u[z]))*scsp.erf(u[z]))
	         um = u[z] - np.sqrt(Zi[z]*np.abs(x)/tau[z])
	         up = u[z] + np.sqrt(Zi[z]*np.abs(x)/tau[z])
	         qin = qin + Ni[z]*np.sqrt(Ti[z]*A[0]/A[z])*Ti[z]*1/8*(2/np.sqrt(np.pi)*((5+2*u[z]**2-(3+2*u[z]**2)/u[z]*np.sqrt(np.abs(x)*Zi[z]/tau[z]))*np.exp(-up[z]**2)+(5+2*u[z]**2+(3+2*u[z]**2)/u[z]*np.sqrt(np.abs(x)*Zi[z]/tau[z]))*np.exp(-um[z]**2))+(3/u[z]+12*u[z]+4*u[z]**3+2*x*Zi[z]/tau[z]*(1/u[z]+2*u[z]))*(scsp.erf(up[z])+scsp.erf(um[z])))
	      
	      #Electron contribution
	      q_currp = qip + Ne*Te*np.abs(Iep)*(2+x) 
	      q_currn = qin + Ne*Te*np.abs(Ien)/(1-x)*(2-x)
	      
	      #energy of recombination on dust surface
	      q_currp = q_currp + E_rec/T0*np.abs(Iip)*Ni[0]
	      q_currn = q_currn + E_rec/T0*np.abs(Iin)*Ni[0]
	      q_currp = q_currp*n0*cs0*T0/4
	      q_currn = q_currn*n0*cs0*T0/4
	      q_seep = 3*Wf/T0*Ne*Te*np.abs(dSEEp*Iep)*n0*cs0*T0/4 #Smirnov PPCF 2007
	      zeta = - x*Te*T0/Wf
	      q_seen = 3*Wf/T0*Ne*Te*np.abs(dSEEn*Ien)*n0*cs0*T0/4*(1+zeta)*(1+2*zeta)/(1+3*zeta) #Vignitchouk PPCF 2014
                
	      #Backscatterings and recombination
	      fcr = 0.51*(1-1.74*(Te*T0)**(-0.23)*np.exp(-0.0043*Te*T0))
	      RN = 0.3595*np.tanh(Ti*T0/3.1132) + 0.0246*np.exp(-((Ti*T0-20.8472)/42.1102)**2)
	      RE = 0.4988*np.tanh(Ti*T0/2.0444) + 0.1274*np.exp(-((Ti*T0+3.2303)/46.0212)**2)
	      qEBS = sw_BS*2*np.abs(qe*Te*T0*cs0*Ne*n0*Ie0/4)*np.exp(-x)*fch
	      qIBSp = sw_BS*RE*qip
	      qIBSn = sw_BS*RE*qin
	      #Techniquement à mettre dans l'équation de la chaleur mais la correction est faible (qq %)
	      qRECp = sw_REC*(1-RN)*Iip*cs0*Ne*n0/4*((E_rec+E_mol)*qe-kb*T) 
	      qRECn = sw_REC*(1-RN)*Iin*cs0*Ne*n0/4*((E_rec+E_mol)*qe-kb*T) 
                
	      #Totatal heat flux/4pia^2 excluding the terms depending on Td
	      qred = (q_currp-q_seep-qEBS-qIBSp+qRECp)*HeavisideX + (q_currn-q_seen-qEBS-qIBSn+qRECn)*HeavisidemX


	      #solve for charge
	      sw_sign = np.diff(np.sign(eqtot))
	      inds = [i for i in range(5000-1) if sw_sign[i] != 0]
	      x1 = (x[0:-2] + x[1:-1])/2
	      chi = np.mean(x1[inds])
	      #if len(chi)>1: 
	      #   j = (np.abs(chi - Chi[it])).argmin()
	      #   chi = chi[inds[j]]
	      qred = np.mean(qred[inds]) #scint.interp1d(x, qred, chi)
    		
	      #solve for temperature
	      td = scinteg.odeint(eqTemp, T, t_f, args = (ktemp, qred, chi, sw_th, Te*T0, Wf, Pvap_A, Pvap_B, hsub, mI, epsilon_d, Tw, kb, qe, me, h, sigma_SB, ad, epsilon0, m_d, cp, material))
	      td = td[-1]
	      
	      if isinstance(td, (complex)):
	         print('Dust temperature is not real') 
	         dead = -10
	         break

	      if Td[it]<=meltT and td>=meltT:
	         print('The dust is melting...')
	         qmelt = (m_d*mu/4/np.pi)**(2/3)*(rho_d/3)**1/3*hmel*(1-Psimelt[it])
	         td = scinteg.odeint(eqTemp, T, t_f, args = (ktemp, qred - qmelt, chi, sw_th, Te*T0, Wf, Pvap_A, Pvap_B, hsub, mI, epsilon_d, Tw, kb, qe, me, h, sigma_SB, ad, epsilon0, m_d, cp, material))
	         td = td[-1]
	         psimelt = 1
	         if td<=meltT:
	            td = meltT
	            q0 = qred - q_sub(Td[it], Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te*T0, epsilon0, me, m_d, cp, material) - sw_th*q_th(Td[it], chi, Wf, Te, me, kb, h, qe) - q_bb(Td[it], Tw, epsilon_d, sigma_SB)
	            q1 = qred - q_sub(td, Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te*T0, epsilon0, me, m_d, cp, material) - sw_th*q_th(td, chi, Wf, Te, me, kb, h, qe) - q_bb(Td, Tw, epsilon_d, sigma_SB);
	            q2 = q0 - q1
	            psimelt = Psimelt[it] + q2*(1/hmel)*(3/rho_d)**(2/3)*(4*np.pi/mu/m_d)**(1/3)*dt*(rho_s/cs0)
	      elif Td[it]>=meltT and td<=meltT:
	         print('solidification')
	         qmelt = (m_d*mu/4/pi)**(2/3)*(rho_d/3)**1/3*hmel*(Psimelt[it])
	         td = scinteg.odeint(eqTemp, T, t_f, args = (ktemp, qred - qmelt, chi, sw_th, Te*T0, Wf, Pvap_A, Pvap_B, hsub, mI, epsilon_d, Tw, kb, qe, me, h, sigma_SB, ad, epsilon0, m_d, cp, material))
	         td = td[-1] 
	         psimelt = 0
	         if td>=meltT:
	            td = meltT
	            q0 = qred - q_sub(Td[it], Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te*T0, epsilon0, me, m_d, cp, material) - sw_th*q_th(Td[it], chi, Wf, Te, me, kb, h, qe) - q_bb(Td[it], Tw, epsilon_d, sigma_SB)
	            q1 = qred - q_sub(td, Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te*T0, epsilon0, me, m_d, cp, material) - sw_th*q_th(td, chi, Wf, Te, me, kb, h, qe) - q_bb(Td([it], Tw, epsilon_d, sigma_SB))
	            q2 = q0 - q1
	            psimelt = Psimelt[it] + q2*(1/hmel)*(3/rho_d)**(2/3)*(4*pi/mu/m_d)**(1/3)*dt*(rho_s/cs0)
                
	      if psimelt==1 and not hasmelt:
	         print('The dust has finished melting.') 
	         hasmelt = True
                
	      if Td[it]<=boilT and td>=boilT:
	         print('The dust is boiling...')
	         qboil = (m_d*mu/4/np.pi)**(2/3)*(rho_d/3)**1/3*hsub*(1-Psiboil[it])
	         td = scinteg.odeint(eqTemp, T, t_f, args = (ktemp, qred - qboil, chi, sw_th, Te*T0, Wf, Pvap_A, Pvap_B, hsub, mI, epsilon_d, Tw, kb, qe, me, h, sigma_SB, ad, epsilon0, m_d, cp, material))
	         td = td[-1] 
	         psiboil = 1
	         if td<=boilT:
	            td = boilT
	            q0 = qred - q_sub(Td[it], Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te*T0, epsilon0, me, m_d, cp, material) - sw_th*q_th(Td[it], chi, Wf, Te, me, kb, h, qe) - q_bb(Td[it], Tw, epsilon_d, sigma_SB)
	            q1 = qred - q_sub(td, Pvap_A, Pvap_B, hsub, mI, kb, qe, ad, Te*T0, epsilon0, me, m_d, cp, material) - sw_th*q_th(td, chi, Wf, Te, me, kb, h, qe) - q_bb(Td[it], Tw, epsilon_d, sigma_SB)
	            q2 = q0 - q1
	            psiboil = Psiboil[it] + q2*(1/hsub)*(3/rho_d)**(2/3)*(4*np.pi/mu/m_d)**(1/3)*dt*(rho_s/cs0)
                    
                    
	      #solve for dust mass
	      mu = scinteg.odeint(eqMass, Mu[it], t_f, args = (td, k1, Pvap_A, Pvap_B, mI, kb, Te, hsub, qe, ad, epsilon0, me, mu, m_d, material))
	      mu = mu[-1]
	      if mu<=1e-4:
	         print('The dust has completely evaporated')
	         dead = 4
	         break
                
                
	      #solve for dust trajectory
	      fcent1 = [0, 0, - 2*vphi[it]*vR[it]/R[it]]
	      fcent2 = [vphi[it]**2, 0, 0]
	      fg = [0, - kg, 0]
	      
	      #Solving for R coordinate
	      Rs = scinteg.odeint(eqMot, np.array([R[it], vR[it]]), t_f, args = (Ni, mR, k, mu, Ti, chi, lnLambda, Te, ER, kE, A, Zi, Nn, Tn, Mn, mami, qw, fcent1[0], fcent2[0], fg[0]))
	      
	      #Solving for Z coordinate
	      Zs = scinteg.odeint(eqMot, np.array([Z[it], vZ[it]]), t_f, args = (Ni, mZ, k, mu, Ti, chi, lnLambda, Te, EZ, kE, A, Zi, Nn, Tn, Mn, mami, qw, fcent1[1], fcent2[1], fg[1]))
	      
	      #Solving for toroidal coordinate
	      phis = scinteg.odeint(eqMot, np.array([phi[it], vphi[it]]), t_f, args = (Ni, mphi, k, mu, Ti, chi, lnLambda, Te, Ephi, kE, A, Zi, Nn, Tn, Mn, mami, qw, fcent1[2], fcent2[2], fg[2]))
            
            
	      #Saving the values of the running iteration
	      it = it + 1
	      t[it] = t[it-1] + dt
	      R[it] = Rs[-1, 0]
	      vR[it] = Rs[-1, 1]
	      Z[it] = Zs[-1, 0]
	      vZ[it] = Zs[-1, 1]
	      phi[it] = phis[-1, 0]
	      vphi[it] = phis[-1, 1]
	      Chi[it] = chi
	      Qw[it] = qw
	      Rd[it] = ad
	      Tau[it, :] = tau
	      Td[it] = td
	      Mu[it] = mu
	      IPHId[it] = (np.abs(meshphi - phi[it])).argmin()
	      IRd[it] = (np.abs(meshR[:, 0] - R[it])).argmin()
	      ITHd[it] = (np.abs(meshZ[0, :] - Z[it])).argmin()
	      IZd[it] = 0
	      Nd[it, :] = [Ne, Ni[0]]
	      Ted[it] = Te
	      Md[it, :] = [mR, mZ, mphi]
	      Ed[it,:] = [ER, EZ, Ephi]
	      Psimelt[it] = psimelt
	      Psiboil[it] = psiboil
	      Nnd[it] = Nn
	      Tnd[it] = Tn
	      Rho_d[it] = rho_d
	      Cp[it] = cp
	      
	      try:
	         pass
	         #[R, Z, vR, vZ] = wall_reflection() #Rwall/rho_s, Zwall/rho_s, R, Z, vR, vZ, it
	      except NameError:
	         if R[it]>max(meshR[:, 0]) or R[it]<min(meshR[:, 0]) or Z[it]>max(meshZ[0, :]) or Z[it]<min(meshZ[0, :]):
	            print('Dust outside simulation domain')
	            dead = 1
	            break
	         
	      if it>Ntmax:
	         dead = -1
	         break
	         
	   #Saving the dust trajectory in the temp files
	   t = t[:it+1]
	   R = R[:it+1]
	   vR = vR[:it+1]
	   Z = Z[:it+1]
	   vZ = vZ[:it+1]
	   phi = phi[:it+1]
	   vphi = vphi[:it+1]
	   Chi = Chi[:it+1]
	   Qw = Qw[:it+1]
	   Tau = Tau[:it+1, :]
	   Td = Td[:it+1]
	   Mu = Mu[:it+1]
	   IPHId = IPHId[:it+1]
	   IRd = IRd[:it+1]
	   ITHd = ITHd[:it+1]
	   IZd = IZd[:it+1]
	   Nd = Nd[:it+1, :]
	   Ted = Ted[:it+1]
	   Md = Md[:it+1, :]
	   Ed = Ed[:it+1, :]
	   Psimelt = Psimelt[:it+1]
	   Psiboil = Psiboil[:it+1]
	   Nnd = Nnd[:it+1]
	   Tnd = Tnd[:it+1]
	   Rho_d = Rho_d[:it+1]
	   Cp = Cp[:it+1]
	   dead = -4 #A supprimer après
	   fnameout = folder1 + fname0 + '_' + f"{kkk+1:03d}" + '.mat'
	   scio.savemat(fnameout, {'NKU' : NKU, 'n0' : n0, 'T0' : T0, 'B0' : B0, 'cs0' : cs0, 'rho_s' : rho_s, 'a' : a, 'x0' : x0, 'v0' : v0, 'R' : R, 'vR' : vR, 'Z' : Z, 'vZ' : vZ, 'phi' : phi, 'vphi' : vphi, 'Chi' : Chi, 'Qw' : Qw, 'Rd' : Rd,'Tau' : Tau, 'Td' : Td, 'Mu' : Mu, 'ITHd' : ITHd, 'IRd' : IRd, 'IZd' : IZd, 'IPHId' : IPHId, 'Nd' : Nd, 'Ted' : Ted, 'Md' : Md, 'Ed' : Ed, 'material' : material, 't' : t, 'dt' : dt, 'test' : test, 'dead' : dead, 'Psimelt' : Psimelt, 'Psiboil' : Psiboil, 'Nnd' : Nnd, 'Tnd' : Tnd, 'Rho_d' : Rho_d, 'Cp' : Cp})
	   print('The results are saved in '+fnameout)
        
	   #End of computation phase
	   t3 = time.time()
	   times.append(t3)
	   tcomp = t3 - t2
	   print('-------------------------------------------')
	   print('Computation time for solving dust trajectory: '+f"{tcomp:.2f}"+' s')
	   print('-------------------------------------------')
	   print('-------------------------------------------')
  
	   # Launching the clock for the final steps
	   t4 = time.time()
	   fname1 = [fich for fich in os.listdir(folder1) if fich.startswith(fname0)]
    
	   #Saving all the trajectory data in one big file
	   if len(fname1) == kkktot:
	      try:
	         del fnameout_f
	      except NameError:
	         fnameout_f = folder1+'ongoing_'+fname0+'.txt'
	         scio.savemat(fnameout_f, {'fnameout_f' : fnameout_f}) #place holder
	         print('Computing neutral atom source...')
	         NAS = reassamble_dust_results(folder1, folder2, fname0, rd, mI/mi)
	         
	         scio.loadmat(folder2+fname0+'.mat')
	         #output_data = data_HFPS() Ecrire une fonction qui donne le temps de vie des poussières, la position d'évaporation, et la quantité de tungstène relâché
	         os.remove(fnameout_f)
	      else:
	         print('Reassembling of trajectory already ongoing ')
	   else:
	      print('Total number of file not enough for reassembling')
	
    
	   #Displaying total computation time
	   t5 = time.time()
	   ttot = t5 - t0
	   print('-------------------------------------------')
	   print('Total computation time: '+f"{ttot:.2f}"+' s')

	#Visualization of the trajectory
	print('Plotting dust trajectory...')
	plt.figure(figsize=(5, 5))
	fontsize = 8
	
	#image of the electronic temperature background
	plt.imshow(np.squeeze(fieldTe*T0), cmap = 'turbo', origin = 'lower')
	plt.colorbar(cmap = 'turbo', label = r"$T_e$ [eV]")
	
	#Adding dust trajectories
	data = scio.loadmat(folder2+fname0+'.mat')
	R = data['R']
	Z = data['Z']
	nb_ligns_R = np.shape(R)[0]
	nb_rows_R = np.shape(R)[1]
	R_im = np.reshape(np.zeros(nb_ligns_R*nb_rows_R), (nb_ligns_R, nb_rows_R))
	Z_im = np.reshape(np.zeros(nb_ligns_R*nb_rows_R), (nb_ligns_R, nb_rows_R))
	for i in range(np.shape(R)[0]):
	   for j in range(np.shape(R)[1]):
	      try:
	         R_im[i, j] = (np.abs(meshR[:, 0] - R[i, j])).argmin()
	         Z_im[i, j] = (np.abs(meshZ[0, :] - Z[i, j])).argmin()
	         j_max = j
	      except IndexError:
	         break
	   plt.plot(R_im[i, :j_max], Z_im[i, :j_max], color = 'gray', linewidth = 2.5, label = f"Dust number {i+1}")
	
	#Plot the legends and all things to improve it 
	try:
	   plt.plot(Rwall, Zwall, color = 'black', linewidth = 5)
	except NameError:
	   pass
	plt.axis('equal')
	plt.xlabel('R (m)', fontsize = fontsize)
	plt.ylabel('Z (m)', fontsize = fontsize)
	pos = np.array([25, 75, 125, 175, 225, 275])
	plt.xticks(pos, [f"{meshR[i, 0]*rho_s:.2f}" for i in pos], fontsize = fontsize)
	plt.yticks(pos, [f"{meshZ[0, i]*rho_s:.3f}" for i in pos], fontsize = fontsize)
	plt.legend()
	plt.show(block=False)
	
	print('DUMBO flies freely toward the sky !')
	print('*******************************************')
	print('*******************************************')
