function q_sub=q_sub(T,Pvap_A,Pvap_B,hsub,mI,kb,qe)
q_sub=Pvap(T,Pvap_A,Pvap_B)./sqrt(2.*pi.*mI.*kb.*T).*hsub./qe*(mI); %DUMBO
end
