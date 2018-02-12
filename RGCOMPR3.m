function [M,TSTR,MFP] = RGCOMPR3(Tt,f,TSPR)
[ht,Prt,phit,cpt,Rt,Gamma_airt,at] = FAIR1(f,Tt);
Pr = Prt/TSPR;
[h,Pr,phi,cp,R,Gamma_air,a] = FAIR1(f,T);
V_2 = sqrt(2*(ht-h)*gc);
if (V_2 < 0)
    M = 0;
    T = Tt
else
    M = sqrt(V_2)/a;
end
[TSTR,TSPR,MFP] = MASSFP(Tt,f,M)
end