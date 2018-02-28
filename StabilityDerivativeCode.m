%% StabilityDerivatives====================================================
%Variables=================================================================
S_vtpr
S_ref%Referance Area
sweep;%Sweep Angle
z_w
d
AR%Aspect Ratio
TR%Taper Ratio
x_wing
z
Cm_acwing
T
z_t
q
S
tau
Cl_alphaw
Cl_betabasic
Cl_betaD
CL_alpha
Lambda%dihedreal angle
CL_alphaVT
S_VT%Vertical Tail Area
z_v
l_VT
L_VT
N_power
N_wing
N
b
CL
CD
x
V_VT%Vertical Tail Volume
c_r%RootChord
%WingGeometery=============================================================
c_bar       = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));
vwt         = (0.724+((3.06*(S_vtpr/S_ref))/(1+cos(sweep)))+(0.4*(z_w/d))+(0.009*AR))
%Longitudinal==============================================================
%N           = (L*cos(AOA))+(D*sin(AOA));
%C           = (D*cos(AOA))-(L*sin(AOA));
Cm_cg       = (-CL*(x_wing/c_bar))+(CD*(z/c_bar))+Cm_acwing+((T*z_t)/(q*S*tau));
Cm_alpha    = ((-x_wing/c_bar)*Cl_alphaw);
%Lateral===================================================================
Cl_betaLam  = -0.25*CL_alpha*Lambda*((2+(1+2*TR))/(3*(1+TR)));
Cl_betawing = Cl_betabasic+Cl_betaD+Cl_betaLam
Cl_betaVT   = (-CL_alphaVT*(S_VT/S_ref)*(z_v/b)*(vwt));
Cl_beta     = Cl_betawing+Cl_betaVT;
%Directional===============================================================
N           = l_VT*L_VT+N_power+N_wing;
Cn          = N/(q*S_ref*b);
Cn_betawing = (CL^2*((1/4*pi*AR)-(tan(sweep)/pi*AR*(AR+4*cos(sweep)))*(cos(sweep)-(AR/2)-(AR^2/8*cos(sweep))+((6*x/c_bar)*(sin(sweep)/AR)))));
Cn_betaVT   = (V_VT*CL_alphaVT*(vwt));
Cn_beta     = Cn_betawing+Cn_betaVT;