%% LongitudinalStability===================================================
%% Variables===========================================================
x_cg                                                             %CG x-axis
x_acwf                                            %AerodynamicCenter x-axis
x_mg
x_wing
x_ach
z
z_D
z_mg
z_cg
I_yymg
theta_ddot
aoa_h
miu                                         %Ground Friction Coefficient
V_r                                                  %Real Velocity
N
m
W
a
    %% WingGeometry========================================================
c_r                                                             %Root Chord
c_bar = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));  %Nicolai Pg580
S                                                            %Wing Area
S_ref                                                %Reference Area
S_h                                                 %H. Tail Area
V_h
TR                                                          %Taper Ratio
mac                                           %Mean Aerodynamic Chord
    %% AeroCoefficients====================================================
Cl
Cl_aoaw
Cl_0wf
Cl_TO
Cl_h
Cl_aoa
Cl_aoah
Cl_1
Cl_0
D
Cd
Cd_TO
Cm_acw
Cm_acwf
    %% PropulsionCharacteristics===========================================
T
z_t
    %% Atmospheric=========================================================
q
q_bar
rho
dynpres_h
%% StabilityDerivativeCalculations=========================================
xbar_cg    = x_cg/c_bar;
xbar_acwf  = x_acwf/c_bar;
Cm_cg      = (-Cl*(x_wing/c_bar))+(Cd*(z/c_bar))+Cm_acw+((T*z_t)/(q*S*c_bar)); %Nicolai Eq 21.3
Cm_aoa     = ((-x_wing/c_bar)*Cl_aoaw); %Nicolai 21.6
Cm_0       = Cm_acwf+(Cl_0wf*(xbar_cg-xbar_acwf));
%% PitchControlSurfaceSizing===============================================
be_frac    = input('Please input the Elevator to wing span ratio:');
dele_max   = input('Please input the Maximum Elevator Deflection Angle [degrees]: ');
dele_maxrd = dele_max*(pi/180);
L_wf       = 0.5*rho*(V_r^2)*Cl_TO*S_ref;
D_TO       = 0.5*rho*(V_r^2)*Cd_TO*S_ref;
Mac_wf     = 0.5*rho*(V_r^2)*Cm_acwf*S_ref*mac;
a_TO       = (T-D-(miu*N))/m;
M_w        = (W*(x_mg-x_cg));
M_D        = (D*(z_D-z_mg));
M_T        = (T*(z_t-z_mg));
ML_wf      = (L_wf*(x_mg-x_acwf));
M_a        = (m*a*(z_cg-z_mg));
L_h        = ((ML_wf+Mac_wf+M_a-M_w+M_D-M_T-(I_yymg*theta_ddot))/(x_ach-x_mg))
Cl_h       = ((2*L_h)/(rho*(V_r^2)*S_h));
tau_e      = (((Cl_h/Cl_aoah)-aoa_h)/dele_maxrd);
ce_frac    = ((1.1873*(tau_e^2))-(0.1602*tau_e)+0.0678);                     %Ratio Equation Approximated in Excel
Cm_dele    = -Cl_aoah*dynpres_h*V_h*be_frac*tau_e;
Cl_dele    = Cl_aoah*dynpres_h*(S_h/S)*be_frac*tau_e;
Cl_hdele   = Cl_aoah*tau_e;
dele_trim  = -((((((T*z_t)/(q_bar*S*mac))*Cl_aoa)+(Cl_1-Cl_0)*Cm_aoa))/((Cl_aoa*Cm_dele)-(Cm_aoa*Cl_dele)))
%%  Plots====================================================================

dele = -25:5:25;
aoa = -4:2:24; 
    for n = 1:length(dele)
        for j = 1:length(aoa)
            AOA(j) = (aoa(j)*(pi/180));
            Cm(n,j)  = (Cm_0+(Cm_aoa*AOA(j))+(Cm_dele*dele(n)));
        end
    end