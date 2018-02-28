function [TSTR,TSPR,MFP] = MASSFP(Tt,f,M)%TSTR is Total to static temperature ratio
gc = 32.174;
[ht,Prt,Phi,Cpt,Rt,Gamma_air,a] = FAIR1(f,Tt);
V = M*a/(1+((Gamma_air-1)/2)*M^2);
Gate = 1;
iter = 0;
while (Gate==1)
    if iter>50
        error('No solution')
    end
h = ht - V^2/(2*778*gc);
[T,Pr,Phi,cp,R,Gamma_air,a] = FAIR2(f,h);
Vn = M*real(a);
if V ~= 0
    Verror = (V-Vn)/V;
else Verror = V-Vn;
end
if abs(Verror)>0.0001
    V=Vn;
    iter = iter+1;
else
    Gate=0;
end
end
TSTR = Tt/T;
TSPR = Prt/Pr;
MFP = (M/TSPR)*sqrt((Gamma_air*gc/R)*(Tt/T));
end