function dydt=eqMass(t,y,td,k1,Pvap_A,Pvap_B,mI,kb,Te,hsub,qe,ad,epsilon0,me,mu,m_d,material)
gamma0 = Pvap(td,Pvap_A,Pvap_B,material)./sqrt(2.*pi.*mI.*kb.*td).*mI;
gamma_vs = vapor_shielding(Te,hsub,mI,kb,qe,ad,epsilon0,me);

if gamma0<gamma_vs
	dydt=-Pvap(td,Pvap_A,Pvap_B,material)./sqrt(2.*pi.*mI.*kb.*td).*k1.*y^(2/3);
else
	dydt=-gamma_vs/m_d*rho_s/cs0;
end
end

