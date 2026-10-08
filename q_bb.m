function q_bb=q_bb(T,Tw,epsilon_d,sigma_SB,qe)
q_bb=sigma_SB*epsilon_d.*(T.^4-Tw.^4)./qe;% last factor to convert to eV/m^2/s
end