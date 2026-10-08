function q_th=q_th(T,x,Wf,Te,me,kb,h,qe)
%Te in eV
    Ithp=16.*pi.*me*(kb.*T).^2./h.^3.*exp(-Wf.*qe./kb./T); 

if x>0
    q_th=(2.*T.*kb./qe).*Ithp./4;
elseif x<=0
    q_th=(T.*kb./qe).*Ithp./4.*...
        (2-2*x.*Te.*qe./kb./T+(x.*Te.*qe./kb./T).^2).*exp(x.*Te.*qe./kb./T);
end

end