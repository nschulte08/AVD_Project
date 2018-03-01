%% StabilityDerivatives====================================================
%Variables=================================================================
S_vtpr%Vertical Tail area extended to fuselage
S_ref%Referance Area
sweep;%Sweep Angle
z_w%z-axis root chord distance from centerline
d%Minimum fuselage depth
AR%Aspect Ratio
TR%Taper Ratio
x_wing
z
Cm_acwing
T%Thrust
z_t%Thrust vector distance to centerline z-axis
q
S
tau
Cl_alphaw%Wing-body lift curve slope
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
c_bar       = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));  %Nicolai Pg580
vwt         = (0.724+((3.06*(S_vtpr/S_ref))/(1+cos(sweep)))+(0.4*(z_w/d))+(0.009*AR)) %Nicolai 21.15
%Longitudinal==============================================================
%N           = (L*cos(AOA))+(D*sin(AOA));
%C           = (D*cos(AOA))-(L*sin(AOA));
Cm_cg       = (-CL*(x_wing/c_bar))+(CD*(z/c_bar))+Cm_acwing+((T*z_t)/(q*S*c_bar)); %Nicolai Eq 21.3
Cm_aoa    = ((-x_wing/c_bar)*Cl_alphaw); %Nicolai 21.6
%Lateral===================================================================
Cl_betaLam  = -0.25*CL_alpha*Lambda*((2+(1+2*TR))/(3*(1+TR))); %Nicolai 21.12
Cl_betawing = Cl_betabasic+Cl_betaD+Cl_betaLam %Nicolai 21.11
Cl_betaVT   = (-CL_alphaVT*(S_VT/S_ref)*(z_v/b)*(vwt)); %Nicolai 21.14
Cl_beta     = Cl_betawing+Cl_betaVT; %Nicolai 21.10
%Directional===============================================================
N           = l_VT*L_VT+N_power+N_wing; %Nicolai 21.18
Cn          = N/(q*S_ref*b); %Nicolai 21.19
Cn_betawing = (CL^2*((1/4*pi*AR)-(tan(sweep)/pi*AR*(AR+4*cos(sweep)))*(cos(sweep)-(AR/2)-(AR^2/8*cos(sweep))+((6*x/c_bar)*(sin(sweep)/AR))))); %Nicolai 21.22
Cn_betaVT   = (V_VT*CL_alphaVT*(vwt)); %Nicolai 21.21
Cn_beta     = Cn_betawing+Cn_betaVT; %Nicolai 21.20