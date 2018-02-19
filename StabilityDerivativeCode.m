%% StabilityDerivatives====================================================
%WingGeometry==============================================================
TR          =  
c_r         = 
c_bar       = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));
z           = 
z_t         = 
tau         =  
%Longitudinal==============================================================
N           = (L*cos(AOA))+(D*sin(AOA));
C           = (D*cos(AOA))-(L*sin(AOA));
Cm_cg       = (-CL*(x_wing/c_bar))+(CD*(z/c_bar))+Cm_acwing+((T*z_t)/(q*S*tau));
Cm_alpha    = ((-x_wing/c_bar)*Cl_alphaw);
%Lateral===================================================================
Cl_betawing = Cl_betabasic+Cl_betaD+Cl_betaLam
Cl_betaLam  = -0.25*CL_alpha*Lambda*((2+(1+2*TR))/(3*(1+TR)));
Cl_betaVT   = (-CL_alphaVT*(S_VT/S_ref)*(z_v/b)*((1+(dalpha/dbeta))*(q_VT/q)));
Cl_beta     = Cl_betawing+Cl_betaVT;
%Directional===============================================================
N           = l_VT*L_VT+N_power+N_wing;
Cn          = N/(q*S_ref*b);
Cn_betawing = (CL^2*((1/4*pi*AR)-(tan(sweep)/pi*AR*(AR+4*cos(sweep)))*(cos(sweep)-(AR/2)-(AR^2/8*cos(sweep))+((6*x/c_bar)*(sin(sweep)/AR))));
Cn_betaVT   = (V_VT*CL_alphaVT*((1+(dalpha/dbeta))*(q_VT/q)));
Cn_beta     = Cn_betawing+Cn_betaVT;