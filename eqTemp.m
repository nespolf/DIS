function dydt=eqTemp(t,y,ktemp,qred,x,sw_th,Te,Wf,Pvap_A,Pvap_B,hsub,mI,epsilon_d,Tw,kb,qe,me,h,sigma_SB)
%all heat fluxes are devided by 4piad^2, units of eV/m^2/s
dydt=ktemp.*(qred-q_sub(y,Pvap_A,Pvap_B,hsub,mI,kb,qe)-sw_th.*q_th(y,x,Wf,Te,me,kb,h,qe)-q_bb(y,Tw,epsilon_d,sigma_SB,qe));
end






