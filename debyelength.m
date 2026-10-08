function lD=debyelength(Te,ne)
%n0 in m^-3 Te in eV B in T
%rs=1.02.*sqrt(2).*sqrt(T0)./(B0.*1e4);%m
lD = sqrt(8.8542e-12*Te./(1.6022e-19.*ne));
%lD=lD./rs;% in rhos units
end