clear all
mainfolder='/Users/fnespoli/scripts/DIS/';
inputfile='/Users/fnespoli/run_dis/input_W7X_AEB30.txt'

    
    addpath(mainfolder)
    tic;


    FID = fopen(inputfile);
    data = textscan(FID,'%s','delimiter','\n');
    a= string(data{:});
    fclose(FID);
    
    rd=str2num(a(2));
    material=data{1}{4};
    x0=str2num(a(6));
    v0=str2num(a(8));
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

load(fnamepbg)
A=1;
    cs0=9.79e3*sqrt(2*T0/A(1));
    rho_s=1.02*sqrt(T0*A(1))/(B0.*1e4); %m
     x0=[x0(1)/rho_s x0(2)/rho_s x0(3)];% in rho_s
    X0=x0(1).*cos(x0(3));  Y0=x0(1).*sin(x0(3)); Z0=x0(2);
    v0=v0./cs0; %in c_s
    v00=sqrt(v0(1).^2+v0(2).^2+v0(3).^2);
    v01=v0;

    v0x=[v0(1).*cos(atan2(v0(1),v0(3))) v0(1).*sin(atan2(v0(1),v0(3))) v0(2)];
    a=[-1 0 0 ]' ; 
    b=v0x'./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2);
    R=2.*((a+b)*(a+b)')/((a+b)'*(a+b))-eye(3);
figure
    for i=1:50
            alpha=180+10.*randn(1);% normally distributed
            beta=rand(1)*2*180-180; %uniformily distributed

            v0(1,i)=v00.*sind(alpha).*cosd(beta); %vertical case
            v0(2,i)=v00.*cosd(alpha);
            v0(3,i)=v00.*sind(alpha).*sind(beta);    
            X00=X0+5e-3./rho_s.*randn(1);
            Y00=Y0+5e-3./rho_s.*randn(1);
            Z00=Z0;
          %  x0=[sqrt(X00.^2+Y00.^2) Z0 atan(Y00/X00)];
            
            v0(1,i)=v00.*cosd(alpha); %horizontal case
            v0(2,i)=v00.*sind(alpha).*cosd(beta);
            v0(3,i)=v00.*sind(alpha).*sind(beta);     

            X00=X0;
                Y00=Y0+5e-3./rho_s.*randn(1);
                Z00=Z0+5e-3./rho_s.*randn(1);


            vp=sqrt(v0(3,i).^2+v0(1,i).^2);
            v1x=[vp.*cos(atan2(v0(3,i),v0(1,i))) vp.*sin(atan2(v0(3,i),v0(1,i))) v0(2,i)]';
            v1R=R*v1x;
                        plot3([0 v1R(1)],[0 v1R(2)],[0 v1R(3)]); hold on
            xlabel('vx'); zlabel('vy'); ylabel('vz')

            vperp=sqrt(v1R(1).^2+v1R(2).^2);
            v1=[vperp.*cos(atan2(v1R(2),v1R(1))-x0(3)) v1R(3) vperp.*sin(atan2(v1R(2),v1R(1))-x0(3))];
            %plot3([0 v0(1,i)],[0 v0(3,i)],[0 v0(2,i)]); hold on
            %scatter3(X00,Y00,Z00); hold on; axis equal

            %plot3([0 v1(1)],[0 v1(3)],[0 v1(2)]); hold on
            xlabel('vR'); zlabel('vZ'); ylabel('vPhi')

    end

%hold on;             plot3([0 v01(1)],[0 v01(3)],[0 v01(2)],'k','linewidth',3); hold on
hold on;             plot3([0 v0x(1)],[0 v0x(2)],[0 v0x(3)],'k','linewidth',3); hold on
    %% rotation with Euler angles https://en.wikipedia.org/wiki/Rotation_matrix
    v0x=[v01(1).*cos(atan2(v01(1),v01(3))) v01(1).*sin(atan2(v01(1),v01(3))) v01(2)];
    v0x=v0x./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2);
    a=acos(-v0x(1)./sqrt(v0x(1).^2+v0x(2).^2));
    b=acos(v0x(3)./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2));b=0;
    c=acos(v0x(3)./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2)); c=0;

Rot=[cos(b).*cos(c) sin(a).*sin(b).*cos(c)-cos(a).*sin(c) cos(a).*sin(b).*cos(c)+sin(a).*sin(c);
    cos(b).*sin(c) sin(a).*sin(b).*sin(c)+cos(a).*cos(c) cos(a).*sin(b).*sin(c)-sin(a).*cos(c);
    -sin(b) sin(a).*cos(b) cos(a).*cos(b)];
x=[1 0 0]';
y=Rot*x;

figure; plot3([0 v0x(1)],[0 v0x(2)],[0 v0x(3)]);
hold on; plot3([0 y(1)],[0 y(2)],[0 y(3)])
xlabel('x'); ylabel('y'); zlabel('z')
%%
    v0x=[v01(1).*cos(atan2(v01(1),v01(3))) v01(1).*sin(atan2(v01(1),v01(3))) v01(2)];
a=acos(v0x(1)./sqrt(v0x(1).^2+v0x(2).^2));
b=acos(v0x(3)./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2));
b=acos(v0x(3)./sqrt(v0x(1).^2+v0x(3).^2));

c=b;
c=-pi/2;

x=[1 0 0]';
Rz=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
Rx=[1 0 0; 0 cos(b) -sin(b); 0 sin(b) cos(b)];
Ry=[cos(c) 0 sin(c); 0 1 0; -sin(c) 0 cos(c)];

y=Rz*x;
y=Rx*y;
%y=Ry*y;
figure; plot([0 v0x(1)]./sqrt(v0x(1).^2+v0x(2).^2),[0 v0x(2)]./sqrt(v0x(1).^2+v0x(2).^2)); hold on; plot([0 y(1)],[0 y(2)],'--')
figure; plot3([0 v0x(1)]./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2),[0 v0x(2)]./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2),[0 v0x(3)]./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2));
hold on; plot3([0 y(1)],[0 y(2)],[0 y(3)],'--')
xlabel('x'); ylabel('y'); zlabel('z')

%% using rodriguez rotation formula https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    v0x=[v01(1).*cos(atan2(v01(1),v01(3))) v01(1).*sin(atan2(v01(1),v01(3))) v01(2)];

a=[1 0 0 ]' ; 
  b=v0x'./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2);

R=2.*((a+b)*(a+b)')/((a+b)'*(a+b))-eye(3);
y=R*a;
figure;
plot3([0 a(1)],[0 a(2)],[0 a(3)]); hold on
plot3([0 b(1)],[0 b(2)],[0 b(3)])
plot3([0 y(1)],[0 y(2)],[0 y(3)],'--')
%% OKKKKK
clear all
mainfolder='/Users/fnespoli/scripts/DIS/';
inputfile='/Users/fnespoli/run_dis/input_W7X_AEB30.txt'

    
    addpath(mainfolder)
    tic;


    FID = fopen(inputfile);
    data = textscan(FID,'%s','delimiter','\n');
    a= string(data{:});
    fclose(FID);
    
    rd=str2num(a(2));
    material=data{1}{4};
    x0=str2num(a(6));
    v0=str2num(a(8));
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

load(fnamepbg)
A=1;
    cs0=9.79e3*sqrt(2*T0/A(1));
    rho_s=1.02*sqrt(T0*A(1))/(B0.*1e4); %m
     x0=[x0(1)/rho_s x0(2)/rho_s x0(3)];% in rho_s
    X0=x0(1).*cos(x0(3));  Y0=x0(1).*sin(x0(3)); Z0=x0(2);
    v0=v0./cs0; %in c_s
    v00=sqrt(v0(1).^2+v0(2).^2+v0(3).^2);
    v01=v0;

    v0x=[v0(1).*cos(atan2(v0(1),v0(3))) v0(1).*sin(atan2(v0(1),v0(3))) v0(2)];
    a=[-1 0 0 ]' ; 
    b=v0x'./sqrt(v0x(1).^2+v0x(2).^2+v0x(3).^2);
    b=v0'./sqrt(v0(1).^2+v0(2).^2+v0(3).^2);
    R=2.*((a+b)*(a+b)')/((a+b)'*(a+b))-eye(3);
figure
    for i=1:50%0
            alpha=180+10.*randn(1);% normally distributed
            beta=rand(1)*2*180-180; %uniformily distributed

           clear v1
            v1(1,1)=v00.*cosd(alpha); %horizontal case
            v1(2,1)=v00.*sind(alpha).*cosd(beta);
            v1(3,1)=v00.*sind(alpha).*sind(beta);    
           
            v1R=R*v1;
            plot3([0 v1R(1)],[0 v1R(3)],[0 v1R(2)]); hold on
                        xlabel('vR'); zlabel('vPhi'); ylabel('vZ')


    end

hold on;             plot3([0 v01(1)],[0 v01(3)],[0 v01(2)],'k','linewidth',3); hold on
