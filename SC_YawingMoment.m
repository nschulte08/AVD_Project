%% DirectionalStability====================================================
%% Variables===============================================================
x
K_f1
K_f2
U1
    %% WingGeometry========================================================
S
S_v
S_vtpr
S_ref
S_s
sweep
z_w
b
d
l_vt
AR
TR
V_vt
c_r
c_bar       = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));
    %% AeroCoefficients====================================================
Cl_aoavt
Cd_y
    %% Atmospheric=========================================================
rho
dynpres_v
    %% Commonterms=========================================================
vwt         = (0.724+((3.06*(S_vtpr/S_ref))/(1+cos(sweep)))+(0.4*(z_w/d))+(0.009*AR)) %Nicolai 21.15
%% StabilityDerivativeCalculations=========================================
Cn_betawing = (Cl^2*((1/4*pi*AR)-(tan(sweep)/pi*AR*(AR+4*cos(sweep)))*(cos(sweep)-(AR/2)-(AR^2/8*cos(sweep))+((6*x/c_bar)*(sin(sweep)/AR))))); %Nicolai 21.22
Cn_betavt   = (V_vt*Cl_aoavt*(vwt)); %Nicolai 21.21
Cn_beta     = Cn_betawing+Cn_betavt; %Nicolai 21.20
Cn_0        = 
%% YawControlSurfaceSizing=================================================
    %% Asymmetric Thrust===================================================
br_frac     = 1
delr_max    = 30;
V_conmin    = 16 
N_A         = -T_L*y_T;
Cn_delr     = N_A/(q*S*b*delr_max);
tau_r       = (Cn_delr*b_vt/(-CL_aoav*V_v*dynpres_v*br_frac));
cr_frac     = ((1.1873*(tau_r^2))-(0.1602*tau_r)+0.0678);
    %% Cross-wind Landing==================================================
V_v         = ((l_vt*S_v)/(b*S));
V_w         = 28
V_t         = sqrt((U1^2)+(V_w^2));
S_s         = 84
F_w         = 0.5*rho*(V_w^2)*S_s*Cd_y;
br_frac     = 1
cr_frac     = 1
beta        = atan(V_w/U1);
Cn_beta     = K_f1*Cl_aoav*(vwt)*((l_vt*S_v)/S)
Cy_beta     = K_f2*Cl_aoav*(vwt)*((S_v)/S)
tau_r       = ((1.5278*cr_frac^3)-(2.7083*cr_frac^2)+(2.2139*cr_frac)+0.0543);
Cy_delr     = Cl_aoav*dynpres_v*tau_r*br_frac*(S_v/S);
Cn_delr     = -Cl_aoav*V_v*dynpres_v*tau_r*br_frac; 
%% Plots===================================================================
delr = -25:5:25;
beta = -4:2:24; 
    for n = 1:length(dele)
        for j = 1:length(aoa)
            BETA(j) = (beta(j)*(pi/180));
            Cn(n,j) = (Cn_0+(Cn_beta*BETA(j))+(Cn_delr*delr(n))+(Cl_dela*dela));
        end
    end