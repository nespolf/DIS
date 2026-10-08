function [R,Z,vR,vZ]=wall_reflection(Rwall,Zwall,R,Z,vR,vZ,it)
%try
    if ~inpolygon(R(it),Z(it),Rwall,Zwall)
        
        
        
        DR=R(it)-R(it-1); DZ=Z(it)-Z(it-1);
        m=DZ/DR; q=Z(it-1)-m*R(it-1);
        for iw=2:length(Rwall)
            DRw=Rwall(iw)-Rwall(iw-1);
            DZw=Zwall(iw)-Zwall(iw-1);
            m1=DZw/DRw; q1=Zwall(iw-1)-m1*Rwall(iw-1);
            if isinf(m)
                R00=R(it); Z00=m1*R00+q1;
            elseif isinf(m1)
                R00=Rwall(iw); Z00=m*R00+q;
            else
                R00=(q1-q)/(m-m1); Z00=m1*R00+q1;
            end
            if R00<=max(R(it-1:it)) && R00>=min(R(it-1:it)) && Z00<=max(Z(it-1:it)) && Z00>=min(Z(it-1:it))
                indw=[iw-1 iw];  R0=R00; Z0=Z00;
            end
        end
        
        
        
        DRw=Rwall(indw(2))-Rwall(indw(1)); DZw=Zwall(indw(2))-Zwall(indw(1));
        d=sqrt(DR.^2+DZ.^2)-sqrt((R(it)-R0).^2+(Z(it)-Z0).^2);
        alpha1=atan2(DZ,DR);
        alpha2=atan2(DZw,DRw);
        
        alpha3=alpha2-alpha1;
        alpha4=alpha2+alpha3;
        Rf=R0+d*cos(alpha4);
        Zf=Z0+d*sin(alpha4);
        v=sqrt(vR(it).^2+vZ(it).^2);
        vRf=v*cos(alpha4);
        vZf=v*sin(alpha4);
        
        
        R(it)=Rf; Z(it)=Zf; vR(it)=vRf; vZ(it)=vZf;
        disp('powder grain reflected on the wall')
        
        
    end
% catch
%     figure; plot(R,Z)
%     hold on; plot(Rwall,Zwall)
%     axis equal
%     keyboard
% end