%% DirectionalStability====================================================
%% Variables===============================================================
x%distance from the aircraft c.g.to the wing aerodynamic center
K_f1%contribution of the fuselage to the derivative Cn_beta
K_f2%contribution of the fuselage to the derivative Cy_beta 
U1%Forward Velocity
V_conmin%Minimum Control Velocity
V_w%Cross-wind Velocity
    %% WingGeometry========================================================
S%Wing Area
S_v%Vertical Tail Area
S_vtpr%V. Tail Area extended to fuselage
S_ref%Referance Area
S_s%Projected Side Area
sweep%Wing Sweep Angle
z_w%Distance b/w root chord and centerline along z-zxis
b%Wing Span
b_vt%Span of vertical Tail
d%Max. Fuselage Depth
l_vt%x-axis distance b/w cg and V.Tail ac
AR%Aspect Ratio
TR%Taper Ratio
V_vt%Vertical Tail Volume
c_r%Root Chord
c_bar       = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));
    %% AeroCoefficients====================================================
Cl_aoavt%Lift Coefficient due to AOA of the verticcal tail
Cd_y%Side Drag Coefficient
    %% PropusionCharacteristics============================================
T_L%Thrust of the operative engine
y_T%y-axis engine placement from centerline
    %% Atmospheric=========================================================
rho%Air Density 
dynpres_v%Vertical Tail Dynamic Pressure Ratio
    %% Commonterms=========================================================
vwt         = (0.724+((3.06*(S_vtpr/S_ref))/(1+cos(sweep)))+(0.4*(z_w/d))+(0.009*AR)) %Nicolai 21.15
%% StabilityDerivativeCalculations=========================================
Cn_betawing = (Cl^2*((1/4*pi*AR)-(tan(sweep)/pi*AR*(AR+4*cos(sweep)))*(cos(sweep)-(AR/2)-(AR^2/8*cos(sweep))+((6*x/c_bar)*(sin(sweep)/AR))))); %Nicolai 21.22
Cn_betavt   = (V_vt*Cl_aoavt*(vwt)); %Nicolai 21.21
Cn_beta     = Cn_betawing+Cn_betavt; %Nicolai 21.20
Cn_0        = 0 %Just a Placeholder
%% YawControlSurfaceSizing=================================================
    %% Asymmetric Thrust===================================================
br_frac     = input('Please input the Rudder to wing span ratio:');
delr_max    = input('Please input the Maximum Rudder Deflection Angle [degrees]: ');
N_A         = -T_L*y_T;
Cn_delr     = N_A/(q*S*b*delr_max);
tau_r       = (Cn_delr*b_vt/(-CL_aoav*V_v*dynpres_v*br_frac));
cr_frac     = ((1.1873*(tau_r^2))-(0.1602*tau_r)+0.0678);
    %% Cross-wind Landing==================================================
V_v         = ((l_vt*S_v)/(b*S));
V_t         = sqrt((U1^2)+(V_w^2));
F_w         = 0.5*rho*(V_w^2)*S_s*Cd_y;
br_frac1     = input('Please input the Rudder to wing span ratio:');
cr_frac1     = input('Please input the Rudder to chord length ratio:');
beta        = atan(V_w/U1);
Cn_beta     = K_f1*Cl_aoav*(vwt)*((l_vt*S_v)/S)
Cy_beta     = K_f2*Cl_aoav*(vwt)*((S_v)/S)
tau_r       = ((1.5278*cr_frac1^3)-(2.7083*cr_frac1^2)+(2.2139*cr_frac1)+0.0543);
Cy_delr     = Cl_aoav*dynpres_v*tau_r*br_frac1*(S_v/S);
Cn_delr     = -Cl_aoav*V_v*dynpres_v*tau_r*br_frac1; 
%% Plots===================================================================
delr = -25:5:25;
beta1 = -4:2:24; 
    for n = 1:length(dele)
        for j = 1:length(aoa)
            BETA(j) = (beta1(j)*(pi/180));
            Cn(n,j) = (Cn_0+(Cn_beta*BETA(j))+(Cn_delr*delr(n))+(Cl_dela*dela));
        end
    end