v0=5;
for i=1:200
 %   v(i)=random('Normal',v0,1);
    alpha(i)=random('Normal',180,10);
    beta(i)=random('Uniform',-180,180);
end

figure;
hist(alpha);

figure;
hist(beta);

figure; scatter(alpha,beta)

figure; plot(sqrt((sind(alpha).*sind(beta)).^2+(sind(alpha).*cosd(beta)).^2+(cosd(alpha)).^2))

vR=v0.*sind(alpha).*cosd(beta);
vR(2,:)=0.*vR(1,:);
vZ=v0.*cosd(alpha);
vZ(2,:)=0.*vZ(1,:);
vphi=v0.*sind(alpha).*sind(beta);
vphi(2,:)=0.*vphi(1,:);
figure; plot3(vR,vphi,vZ)