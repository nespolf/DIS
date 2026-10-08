clear all 

folder='~/run_DIS/final/';

fname='KSTAR_test_new_B_4e-05m_NAS_EMC3.txt'; material='B';
fname_imp=[folder,'KSTAR_test_new_B_4e-05m_sym_input.imp_',material];

% fname='KSTAR_34002_B_4e-05m_NAS_EMC3.txt'; material='B';
% fname_imp=[folder,'KSTAR_34002_B_4e-05m_sym_input.imp_',material];

%fname='KSTAR_34317_B_4e-05m_NAS_EMC3.txt'; material='B';
%fname_imp=[folder,'KSTAR_34317_B_4e-05m_sym_input.imp_',material];

fname='KSTAR_34317_2e19_B_4e-05m_NAS_EMC3.txt'; material='B';
fname_imp=[folder,'KSTAR_34317_2e19_B_4e-05m_sym_input.imp_',material];

 fname='KSTAR_34317_1e19_B_4e-05m_NAS_EMC3.txt'; material='B';
 fname_imp=[folder,'KSTAR_34317_1e19_B_4e-05m_sym_input.imp_',material];
% 
 fname='KSTAR_34317_3e19_B_4e-05m_NAS_EMC3.txt'; material='B';
 fname_imp=[folder,'KSTAR_34317_3e19_B_4e-05m_sym_input.imp_',material];

 fname='KSTAR_34317_05e19_B_4e-05m_NAS_EMC3.txt'; material='B';
 fname_imp=[folder,'KSTAR_34317_05e19_B_4e-05m_sym_input.imp_',material];
% 
%  fname='KSTAR_34317_15e18_B_4e-05m_NAS_EMC3.txt'; material='B';
%  fname_imp=[folder,'KSTAR_34317_15e18_B_4e-05m_sym_input.imp_',material];

%   fname='KSTAR_34317_25e18_B_4e-05m_NAS_EMC3.txt'; material='B';
%  fname_imp=[folder,'KSTAR_34317_25e18_B_4e-05m_sym_input.imp_',material];

fname='KSTAR_35822_nopuff_W2_B_4e-05m_NAS_EMC3.txt'; material='B';
 fname_imp=[folder,'KSTAR_35822_nopuff_W2_B_4e-05m_sym_input.imp_',material];


a=load([folder,fname]);
sum(a(:,2))
S=a(:,2);
R=a(:,3);
Z=a(:,4);
phi=a(:,5);

R1=linspace(min(R),max(R),floor((max(R)-min(R))./0.5));
dR=mean(diff(R1))/2;
Z1=linspace(min(Z),max(Z),floor((max(Z)-min(Z))./0.5));
dZ=mean(diff(Z1))/2;
for i=1:length(R1)
    for j=1:length(Z1)
        cond=((R>=(R1(i)-dR)).*(R<(R1(i)+dR))).*((Z>=(Z1(j)-dZ)).*(Z<(Z1(j)+dZ)));
        S1(i,j)=sum(S(cond>0));
        R2(i,j)=R1(i);
        Z2(i,j)=Z1(j);
    end
end
S1=S1(:); R2=R2(:); Z2=Z2(:);
S3=S1(S1>0);
R3=R2(S1>0);
Z3=Z2(S1>0);
%%
figure; 
subplot(2,1,2)
scatter(R3,Z3,20,S3,'fill'); axis equal; xlabel('R [m]'); ylabel('Z [m]')
set(gca,'Fontsize',14); grid on; box on; ylim([70 90]) ;xlim([150 190])
subplot(2,1,1)
scatter3(R.*cosd(phi),R.*sind(phi),Z,20,S,'fill','MarkerFaceAlpha',0.2); axis equal; xlabel('x [m]'); ylabel('y [m]');  zlabel('z [m]')
set(gca,'Fontsize',14); grid on; box on;
%
%% write new source
phi1=linspace(-180,180,360/4);
R4=[]; Z4=[]; phi4=[]; S4=[];
for i=1:length(phi1)
R4=[R4 R3'];
Z4=[Z4 Z3'];
S4=[S4 S3'];
phi4=[phi4 phi1(i).*ones(size(S3'))];
end

figure;
scatter3(R4.*cosd(phi4),R4.*sind(phi4),Z4,20,S4,'fill','MarkerFaceAlpha',0.2); axis equal; xlabel('x [m]'); ylabel('y [m]');  zlabel('z [m]')
S4=S4.*sum(S3)./sum(S4);
clear b
b(:,1)=a(1,1).*ones(size(S4));
b(:,2)=S4;
b(:,3)=R4;
b(:,4)=Z4;
b(:,5)=phi4;
b(:,6)=0.*ones(size(S4));
b(:,7)=0.*ones(size(S4));
b(:,8)=0.*ones(size(S4));


write_EMC3_imp_input(b,fname_imp,material)


%% for gas puff
phi1=linspace(-180,180,360/4);
clear B
A=[0.356 0.01 162.0 -120.6 0.0 0 0 0];
for i=1:length(phi1)
    B(i,:)=A; B(i,2)=A(2)/length(phi1); B(i,5)=phi1(i);
end
write_EMC3_imp_input(B,'input.imp_Ne_sym','Ne')
