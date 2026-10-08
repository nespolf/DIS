clear all
%load('/Users/fnespoli/run_DIS/final/DIIID_USN_test_B_7.5e-05m.mat')
load('/Users/fnespoli/run_DIS/final/LHD_160035_B_test_sSB_B_7.5e-05m.mat')
%load('/Users/fnespoli/dust_results_9/final/C_T040_n05e+18_rd5_NKU0_theta30_t3.mat')
A=1; %H plasma
mami=1; %m_atom/m_i
epsilon0=8.8542e-12;
mi=A(1).*1.67e-27;%kg, H
me=9.1094e-31; %kg
Tw=300;
qe=1.6022e-19;
sigma_SB=5.67e-8; %W/m^2/sec/K^4
h=6.6261e-34; %J*s Planck constant
kb=1.3807e-23; %J/K Boltzmann constant
mime=sqrt(mi/me);
cs0=9.79e3*sqrt(T0/A(1));
rho_s=1.02*sqrt(2*T0*A(1))/(B0.*1e4); %mmi=2*1.67e-27;%kg
me=9.1094e-31; %kg
h=6.6261e-34; %J*s Planck constant
kb=1.3807e-23; %J/K Boltzmann constant

%T0=50; B0=1;
cs0=9.79e3*sqrt(T0);

rho_s=1.02*sqrt(2*T0)/(B0.*1e4); %m
rho_d=2267; %graphite kg/m^3

m_d=4*pi/3*(a*rho_s).^3.*rho_d;
k=(3/(4*pi*rho_d))^(2/3)*mi*m_d^(-1/3)*n0*rho_s;
kg=rho_s/cs0.^2*9.81;
kE=(3/4/pi/rho_d/m_d.^2).^(1/3).*(T0/cs0).^2.*(4*pi*epsilon0);
kL=(3/rho_d/m_d.^2*(4*pi).^2).^(1/3).*(T0.*B0.*rho_s/cs0).*(epsilon0);

v=sqrt(vphi.^2+vR.^2+vZ.^2);
M=sqrt(Md(:,:,1).^2 + Md(:,:,2).^2+Md(:,:,3).^2);
u=sqrt((Md(:,:,1)-vR).^2+(Md(:,:,2)-vZ).^2+(Md(:,:,3)-vphi).^2)./sqrt(Ted.*Tau);
u_tilde=u-mean(u,1);
E=sqrt(Ed(:,:,1).^2+Ed(:,:,2).^2+Ed(:,:,3).^2);


u(:,:,1)=(Md(:,:,1)-vR)./sqrt(Ted.*Tau);
u(:,:,2)=(Md(:,:,2)-vZ)./sqrt(Ted.*Tau);
u(:,:,3)=(Md(:,:,3)-vphi)./sqrt(Ted.*Tau);


lD=debyelength(Ted.*T0,Nd.*n0);
ad=(3/4/pi*m_d*Mu/rho_d).^(1/3);

clear Fd Fc Zeta
try
    Tnd=Tnd'; Nnd=Nnd';
end
for i=1:size(Chi,1)
    for j=1:size(Chi,2)
        for l=1:3
            lnLambda=Coul_log(Chi(i,j),Tau(i,j),ad(i,j),u(i,j,l),lD(i,j));
            Fd(i,j,l)=zeta_drag(abs(u(i,j,l)),Chi(i,j),Tau(i,j),lnLambda).*Mu(i,j).^(2/3).*k.*Nd(i,j).*u(i,j,l).*sqrt(Ted(i,j).*Tau(i,j));
        end
        %Fd(i,j)=zeta_drag(abs(u(i,j)),Chi(i,j),Tau(i,j),lnLambda(i,j)).*Mu(i,j).^(2/3).*k.*Nd(i,j).*u(i,j).*sqrt(Ted(i,j).*Tau(i,j));
        Mn=0;
        if exist('Nnd')
            for l=1:3
                Fn(i,j,l)=zeta_drag_n(abs(Mn-u(i,j,l))./sqrt(Tnd(i,j)./mami)).*pi.*Mu(i,j).^(2/3).*k.*mami.*Nnd(i,j).*(Mn-u(i,j,l)).*sqrt(Tnd(i,j)./mami);
            end
        end
        %end
    end
end
Fc(:,:,1)=-(vphi.^2)./R.*Mu;
Fc(:,:,2)=Fc(:,:,1).*0;
Fc(:,:,3)=(vphi.*vR)./R.*Mu;
t=(1:length(R)).*rho_s./cs0.*dt;

Fg=kg.*Mu;
for l=1:3
    FE(:,:,l)=(Chi.*Ted.*Mu.^(1/3).*Ed(:,:,l).*kE);
end


figure;
plot(t,squeeze(Fd(1,:,:)));
hold on
plot(t,squeeze(Fc(1,:,:)));
plot(t,Fg(1,:));
plot(t,squeeze(FE(1,:,:)));
plot(t,squeeze(Fn(1,:,:)));


figure;
plot(t,squeeze(sqrt(Fd(1,:,1).^2+Fd(1,:,2).^2+Fd(1,:,3).^2)));
hold on
plot(t,squeeze(sqrt(Fc(1,:,1).^2+Fc(1,:,2).^2+Fc(1,:,3).^2)));
plot(t,squeeze(sqrt(FE(1,:,1).^2+FE(1,:,2).^2+FE(1,:,3).^2)));
plot(t,Fg(1,:));
plot(t,squeeze(sqrt(Fn(1,:,1).^2+Fn(1,:,2).^2+Fn(1,:,3).^2)));
set(gca,'yscale','log'); grid on; box on



