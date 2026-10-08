function [rho_d cp]= cp_func(material,td,meltT)
if strcmp(material,'B') %https://webbook.nist.gov/
    mmol=10.81;
    alpha_v=3.*6.4e-6; rho_d = 2340; %T. Lundstrom Journal of Alloys and Compunds 267 1998
    if td>=meltT
        rho_d = 2080;
        CP=[31.75 0 0 0 0];
        alpha_v=0;
       % rho_d=(2.468-5.31e-5.*td).*1000; %F. Millot International Journal of Thermophysics, Vol. 23, No. 5, September 2002 ( 2002)
    elseif  td >= 1800
        CP=[ 25.12664 1.975493 0.338395 -0.040032  -2.635578];
    else
        CP=[ 10.18574 29.24415 -18.02137 4.212326 -0.550999];
    end
  
elseif strcmp(material,'B4C') %https://webbook.nist.gov/
    mmol=55.255; rho_d=2520; %solid, cannot find liquid state
    alpha_v=3.*5.65e-6;% G.VTsagareishvili  Journal of the Less Common Metals 117 1986
    rho_d=rho_d./(1+alpha_v.*(td-298));
    if td>=meltT
        CP=[136 0 0 0 0];
    else
        CP=[95.99853	23.16513	-0.409604	0.081414	-4.395208];
    end
    
elseif strcmp(material,'Li') %https://webbook.nist.gov/
    mmol=6.94; 
    alpha_v=3.*46.8e-6;%THERMAL EXPANSION OF LITHIUM, 77° TO 300° K.  WB Pearson · 1954
    alpha_v=1.968e-4-6.724./td.^2+1.413e5./td.^4;% an evaluation of some thermodynamic and transport properties of solid and liquid lithium over the tempreature rnage 200-1700K R.K. Williams 1988
    if td>=meltT
        rho_d=512; 
        alpha_v=1.01e-4./(0.5584-1.01e-4.*td);
        CP=[29 0 0 0 0];
    else
        rho_d=534;
        CP=[169.5520	-882.7110	1977.438	-1487.312	-1.609635];
    end
    
elseif strcmp(material,'BN') %https://webbook.nist.gov/
    mmol=24.82;
    rho_d=2100; %hexagonal, 3450 for c-BN
    rho_d=2250; %from Momentive specs
    alpha_par=38e-6; %hexagonal, M.E. Levinshtein, S.L. Rumyantsev, M.S. Shur Properties of Advanced Semiconductor Materials, 2001
    alpha_perp=-2.7e-6;
    alpha_v=alpha_par+2.*alpha_perp;
    
    if td<1100
        CP=[-1.626426	101.9817		-81.52482		25.65093		-0.223052	];
    else 
        CP=[78.14750	-23.33367	6.844396	-0.699753	-17.04034];
    end
    
elseif strcmp(material,'C') %Butland JNM 1973
    mmol=12;
   alpha_z=27e-6+3.05e-9.*(td-273);%thermal expansion coefficient
   alpha_b=-1.5e-6.*(td<=523) + (1e-8.*td-6.73e-6).*(td<=673).*(td>523)+ (0.45e-8.*td-3.02e-6).*(td<=873).*(td>673)...
   +(0.9e-6).*(td<=1073).*(td>873)+(2.5e-10.*td+0.63e-6).*(td<=1273).*(td>1073)+(0.95e-6).*(td>1273);
  alpha_v=alpha_z+2.*alpha_b;
    rho_d=2267; %graphite kg/m^3

   CP=[0.539 9.11e-6 -90.27 -43449 1.59e7 -1.4437e9];%Butland JNM 1973 cal/g/K
   CP=CP*4.1868*1000;%J/K/kg
        
    
    
else
    disp(['temperature varying Cp for material ',material,' not implemented yet' ])
end

if strcmp(material,'C')
    cpT=td;
   cp=(CP(1)+ CP(2).*cpT + CP(3).*(cpT).^-1 + CP(4).*(cpT).^-2 + CP(5).*(cpT).^-3 +CP(6).*(cpT).^-4); % function for boron

else
cpT = td/1000;
cp=(CP(1)+ (CP(2).*cpT) + (CP(3).*(cpT).^2) + (CP(4).*(cpT).^3) + (CP(5)./(cpT.^2)))./(mmol.*1e-3); % function for boron
end

  rho_d=rho_d./(1+alpha_v.*(td-298));
end