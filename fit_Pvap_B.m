%% fit B vapor pressure
%For vapor pressure fitted from [American Insitute of Physics handbook 3rd ed. p.4-298] 
T=[1620 1750 1900 2076 2291 2550 2870 3275 3802]+273;
P=[1e-3 1e-2 0.1 1  10  100  1e3 10e3 100e3];
xData=T;
yData=log10(P);
ft = fittype(' a./(x)+b; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=abs(xData);
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
x=linspace(300,5000,200);
yFit=feval(fitres,x);
yFit=10.^(yFit);
figure;
scatter(xData,10.^(yData));
hold on; 
plot(x,yFit);
set(gca,'yscale','log'); grid on

%log (P/Pa) = 13.255 - 21370 / (T/K) wikipedia for solid
%Pvap_A=-40181; Pvap_B=14.8;
plot(x,10.^(14.8-40181./x))

%% B4C
%fit B4C vapor pressure
%data form     [H.E. Robson and P.W. Gilles J.Phys.Chem 1964]

T=[2303 2302 2308 2306 2300 2299 2296 2294 2297 2300 2299 2298 2392 2389 2390 2395 2392 2391 2392 2392 2396 2396 2395 2308 2199 2376 2421 2464 2522 2270 2343 2405 2483 2283 2236 2184 2261];
P=1e-6.*101325.*[6.31 4.64 4.65 4.37 4.28 5.88 5.86 5.25 5.76 6.31 6.57 6.35 15.52 15.82 18.20 15.63 14.42 16.11 15.26 14.45 14.10 15.29 14.56 4.79 0.34 12.67 22.02 39.18 80.55 1.75 4.07 12.05 21.30 2.90 1.80 0.900 2.51];
xData=T;
yData=log10(P);
ft = fittype(' a./(x)+b; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=abs(xData);
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
x=linspace(300,5000,200);
yFit=feval(fitres,x);
yFit=10.^(yFit);
%figure;
scatter(xData,10.^(yData));
hold on; 
plot(x,yFit);
set(gca,'yscale','log'); grid on

%% BN
% data from [Lorenz R, Woolcock J (1928) Z Anorg Chemie 176: 283-304]
T=[1695 1745 1740 1800 1820 1890 1895 1930 1970 2045]+273; 
P=[23 26 27 31 43 67 70 79 97 158].*133.322; %conversion from mmHg to Pa
% data from [D. L. Hildenbrand and W. F. Hall J. Phys. Chem. 1963, 67, 4,
% 888?893]
 T=[1945 1971 1978 2023 2042 2045 1998 1987 1926 1912  1860 1925 1978 2022 2043 1995 1957 1920 1893 1941 1984 2022 2043 2076 2162 2114 2056 2001 1966 ];
% P=1e-5.*101325.*[1.01 1.47 1.64 3.44 4.16 4.43 2.37 1.84 0.82 0.55 0.57 1.42 3.10 6.48 8.60 4.47 2.37 1.26 0.78 3.5 7.1 12.6 15.9 24.1 74.0 42.1 19.1 8.2 4.7];
 P=1e-5.*101325.*[4.44 6.47 7.22 15.1 18.3 19.5 10.4 8.10 3.61 2.43 1.27 3.17 6.92 14.4 19.2 9.96 5.28 2.81 1.74 4.30 8.74 15.5 19.6 29.6 91.0 51.8 23.5 10.1 5.78];
% T=[T T1];
% P=[P P1];
xData=T;
yData=log10(P);
ft = fittype(' a./(x)+b; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=abs(xData);
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
x=linspace(300,5000,200);
yFit=feval(fitres,x);
yFit=10.^(yFit);
%figure;
scatter(xData,10.^(yData));
hold on; 
plot(x,yFit);
set(gca,'yscale','log'); grid on
%%
Tboil=[3246 3770 3925 4200];
Pboil=[1 1 1 1].*101325; %Pa,=1atm
hold on; scatter(Tboil,Pboil)

%%
plot(x,10.^(13.255-21370./x))

%% Li
%For vapor pressure fitted from [American Insitute of Physics handbook 3rd ed. p.4-298] 
T=[344 397 456 531 619 730 871 1068 1324 ]+273;
P=[1e-3 1e-2 0.1 1  10  100  1e3 10e3 100e3];
xData=T;
yData=log10(P);
ft = fittype(' a./(x)+b; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=abs(xData);
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
x=linspace(300,5000,200);
yFit=feval(fitres,x);
yFit=10.^(yFit);
figure;
scatter(xData,10.^(yData));
hold on; 
plot(x,yFit);
set(gca,'yscale','log'); grid on

%log (P/Pa) = 13.255 - 21370 / (T/K) wikipedia for solid
%Pvap_A=-40181; Pvap_B=14.8;
plot(x,10.^(14.8-40181./x))

%% J. Am. Chem. Soc. 1918, 40, 1, 45?49
% Hg data form     %https://physics.nyu.edu/kentlab/How_to/ChemicalInfo/VaporPressure/morepressure.pdf
T=[42.0 76.6 120.0 175.6 250.3 355.9]+273;
P=[ 1 10  100  1e3 10e3 100e3]

xData=T;
yData=log10(P);
ft = fittype(' a./(x)+b; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=abs(xData);
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
x=linspace(300,5000,200);
yFit=feval(fitres,x);
%yFit=10.^(yFit);
figure;
scatter(xData,(yData));
hold on; 
plot(x,yFit);
%%
figure
boilT_Hg=629.88;
%C 
T=linspace(300, 5000,1000);

PvapC_A=-40181; PvapC_B=14.8;  boilTC=3925;
logpC=PvapC_A./T +PvapC_B;
plot(T,logpC,'o'); hold on
% a_C=boilTC/boilT_Hg;
% plot(fitres.a.*a_C./T+fitres.b+log10(a_C))
% B
PvapB_A=-28238; PvapB_B=11.98;  boilTB=4200;
logpB=PvapB_A./T +PvapB_B;
plot(T,logpB); hold on
a_C=T/4176;
logpB1=PvapC_A.*a_C./T +PvapC_B+log10(a_C);
plot(T,logpB1); hold on


ylim([-2 10])

%%
Hvap=508; %kJ/mol
 boilT=4200;
 R=8.31%J/K/mol
 P1=1e-5*1e-3;
 boilT1=1./(1/boilT-R.*log(P1/1)/(Hvap*1000));
 
lnpv=Hvap*1000/R.*(1/boilT1-1./T);

hold on; 
plot(T,10.^(lnpv))

%%  temp            density         heatCap         emis            young           poiss           yielS           ultiS           vapPres  
a=[      0.0000000e+00   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   0.0000000e+00
   2.0000000e+02   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08  7.8523563e-187
   3.0000000e+02   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08  7.3001760e-120
   4.0000000e+02   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.6000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   2.2258710e-86
   5.0000000e+02   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   2.7415742e-66
   6.0000000e+02   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   6.7868252e-53
   8.0000000e+02   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   3.7475721e-36
   1.0000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   4.1591061e-26
   1.2000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   2.0693471e-19
   1.5000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   1.0295956e-12
   2.0000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   5.1227127e-06
   2.5000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   5.3407223e-02
   3.0000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   2.5487858e+01
   3.5000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   2.0879221e+03
   4.0000000e+03   2.2500000e+03   1.7090000e+03   7.5000000e-01   7.0000000e+10   1.0000000e-01   1.1000000e+08   6.0000000e+08   5.6852557e+04];
plot(a(:,1),a(:,9));