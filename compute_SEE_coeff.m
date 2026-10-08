%% Boron
%from CRC Handbook of Chemistry and Physics p. 12-125
delta=[1 1.2 1];
E=[50 150 600]; %[eV]
%Graphite
% delta=[1 1 1];
% E=[300 300 300];
r=linspace(1,10,200);
k=1./(1-r./(exp(r)-1));
figure; plot(r,k)
x=E(1);
YD=max(delta)./(1-exp(-r)).*(x./E(2)).^(1-k).*(1-exp(-r.*(x./E(2)).^k)); %Young-Dekker Formula
kk1=dsearchn(YD',1);
x=E(3);
YD=max(delta)./(1-exp(-r)).*(x./E(2)).^(1-k).*(1-exp(-r.*(x./E(2)).^k));
kk2=dsearchn(YD',1);


x=linspace(0,2.*max(E),1000);
YD=max(delta)./(1-exp(-r(kk1))).*(x./E(2)).^(1-k(kk1)).*(1-exp(-r(kk1).*(x./E(2)).^k(kk1))).*(x<=E(2))+...
    max(delta)./(1-exp(-r(kk2))).*(x./E(2)).^(1-k(kk2)).*(1-exp(-r(kk2).*(x./E(2)).^k(kk2))).*(x>E(2));

figure; plot(x,YD);
hold on; scatter(E,delta);
set(gca,'yscale','log');
set(gca,'xscale','log');
ylim([1e-4 1]); xlim([1e-1 1e3])

% fit DUMBO style
ft = fittype(' a+b.*x+c.*x.^2+d.*x.^3+e.*x.^4+f.*x.^5; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
xData=log10(x);
yData=log10(YD);
xData=xData(isnan(yData)<1);
yData=yData(isnan(yData)<1);
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=1./abs(xData).^2;
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
yFit=feval(fitres,xData);
plot(10.^xData,10.^yData,'or'); hold on
plot(10.^xData,10.^yFit,'k','linewidth',2,'HandleVisibility','off')


%% B4C data from [Zameroski IEEE transactions on plasma science 2006]
x=50:50:1000;
YD=[0.3 0.45 0.52 0.54 0.55 0.53 NaN 0.51 NaN 0.48 NaN 0.45 NaN 0.43 NaN 0.41 NaN 0.4 NaN 0.38]; 
x=[0.1 x];
YD=[1e-6 YD];
ft = fittype(' a+b.*x+c.*x.^2+d.*x.^3+e.*x.^4+f.*x.^5; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
xData=log10(x);
yData=log10(YD);
xData=xData(isnan(yData)<1);
yData=yData(isnan(yData)<1);
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
%opts.Weights=1./abs(xData).^2;
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
xFit=linspace(0,3,1000);
yFit=feval(fitres,xFit);
figure
plot(10.^xData,10.^yData,'or'); hold on
plot(10.^xFit,10.^yFit,'k','linewidth',2,'HandleVisibility','off')
set(gca,'yscale','log');
set(gca,'xscale','log');
ylim([1e-4 1]); xlim([1e-1 1e3])


%% BN
%data from [Smirnov JAP 2003];
T=0:10:100;
T(1)=0.9;
YD=0.173.*T.^0.5;
% data from Prokofiev 2010 10.1109/IVESC.2010.5644458 (Pyrolitic BN!)
T1=[125 125*1.5 250 380 500 625 750 750+125/2 875 1000]
YD1=[2.25 2.51 2.75 2.49  2.25 2.11 2 1.9 1.8 1.7];
x=[T T1];
YD=[YD YD1];

ft = fittype(' a+b.*x+c.*x.^2+d.*x.^3+e.*x.^4+f.*x.^5; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
xData=log10(x);
yData=log10(YD);
xData=xData(isnan(yData)<1);
yData=yData(isnan(yData)<1);
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
opts.Weights=1./abs(xData).^2;
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
xFit=linspace(0,3,1000);
yFit=feval(fitres,xFit);
figure
plot(10.^xData,10.^yData,'or'); hold on
plot(10.^xFit,10.^yFit,'k','linewidth',2,'HandleVisibility','off')
set(gca,'yscale','log');
set(gca,'xscale','log');
ylim([1e-4 1]); xlim([1e-1 1e3])
%%
%% Lithium
% data from CRC handbook
T=[85]; YD=0.5; 
% data from [Capace et al APL 2016]
T1=[20 50 105 155 210 250 310 410 505 600 1000]; %last point is fake to keep the curve for going up for T>600
YD1=[0.35 0.5 0.55 0.52 0.5 0.47 0.45 0.4 0.39 0.35 0.25];
x=[T T1];
YD=[YD YD1];

ft = fittype(' a+b.*x+c.*x.^2+d.*x.^3+e.*x.^4+f.*x.^5; ','independent', {'x'}, 'dependent', 'y' );

clear opts
opts = fitoptions( ft );
xData=log10(x);
yData=log10(YD);
xData=xData(isnan(yData)<1);
yData=yData(isnan(yData)<1);
%opts.Lower =        [0   0  ];
%opts.StartPoint =   [(max(xData)./max(yData)) 2];
% %opts.Upper =        [max(yData)    1000   100 ];
opts.Weights=1./abs(xData).^2;
[fitres,gof,output] = fit([xData'],yData', ft,opts);
fitres
xFit=linspace(0,3,1000);
yFit=feval(fitres,xFit);
figure
plot(10.^xData,10.^yData,'or'); hold on
plot(10.^xFit,10.^yFit,'k','linewidth',2,'HandleVisibility','off')
set(gca,'yscale','log');
set(gca,'xscale','log');
ylim([1e-4 1]); xlim([1e-1 1e3])
