function dydt=eqMass(t,y,td,k1,Pvap_A,Pvap_B,mI,kb)

dydt=-Pvap(td,Pvap_A,Pvap_B)./sqrt(2.*pi.*mI.*kb.*td).*k1.*y^(2/3); 

end

