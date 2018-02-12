function [TSTR,TSPR,M] = RGCOMPR5(Tt,f,MFP)%M>1
M = .5
dM = 0.1
[TSTR,TSPR,MFP0] = MASSFP(Tt,f,M)
Gate = 1;
while Gate==1
    M = M+dM;
    [TSTR,TSPR,MFPn] = MASSFP(Tt,f,M)
    MFPerror = abs(MFPn-MFP0);
    if MFPerror>0.00001
        dM = ((MFP-MFPn)/(MFPn-MFP0))*dM;
        MFP0 = MFPn;
    else Gate=0;
    end
end
end