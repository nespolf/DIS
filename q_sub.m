function q_sub=q_sub(T,Pvap_A,Pvap_B,hsub,mI,kb,qe,ad,Te,epsilon0,me,m_d,cp,material)

q_sub=Pvap(T,Pvap_A,Pvap_B,material)./sqrt(2.*pi.*mI.*kb.*T).*hsub./qe*(mI); %DUMBO

if q_sub>vapor_shielding(Te,hsub,mI,kb,qe,ad,epsilon0,me)*hsub/qe
	q_sub = vapor_shielding(Te,hsub,mI,kb,qe,ad,epsilon0,me)*hsub/qe/m_d/cp;
end

end
