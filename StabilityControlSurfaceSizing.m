%% Control Surface Sizing==================================================
    %By: Christian Allen
    %12/19/18
%% Flight Phases===========================================================
%Class()
%Phase(B)
%% Aircraft Geometry=======================================================
S        = 
TR       = 
c        = 
c_root   = 
b        = 
b_vt     = 
%% Aileron==================================================================
ba_frac   = 
ca_frac   = 
tau_a    = 
fun1     = @(y)(1+(2*((TR-1)/b))*y)*y;
in1      = integral(fun1,y_i,y_o);
Cl_dela  = ((2*CL_aoaw*tau_a*C_root)/(S*b))*(in1);
dela_max = +-25;
Cl       = Cl_dela*dela_max;
L_A      = 2*delL*y_a;
P_ss     = sqrt((2*L_A)/(rho*(S_w+S_vt)*C_Dr*(y_D^3)));
phi1     = ((I_xx)/(rho*(y_D^3)*(S_w+S_vt)*C_Dr))*(log(P_SS^2));
P_dot    = (P_ss^2)/(2*phi1);
%if phi1 > phi_req
t2       = sqrt((2*phi_des)/(P_dot));
%if phi1 < phi_req
t2       = sqrt((2*phi1)/(P_dot))+((phi2*phi1)/(P_ss));
%% Elevators================================================================
be_frac  = 
ce_frac  = 

%% Rudder==================================================================
%Asymmetric Thrust=========================================================
br_frac  = 
delr_max = +-30;
V_conmin = 
N_A      = -T_L*y_T;
Cn_delr  = N_A/(q*S*b*delr_max);
tau_r    = (Cn_delr*b_v/(-CL_aoav*V_v*dynpres_v*b_r));
cr_frac  = %From Fig12.12 Sadrey)
%Cross-wind Landing========================================================
V_w      = 
V_t      = sqrt((U1^2)+(V_w^2));
S_s      = 
F_w      = 0.5*rho*(V_w^2)*S_s*CD_y;
b_frac   = 
c_frac   = 
beta     = atan(V_w/U1);
Cn_beta  = K_f1*CL_aoav*((1-(dalpha/dbeta))*(q_VT/q))*((l_vt*S_v)/S)
Cy_beta  = K_f2*CL_aoav*((1-(dalpha/dbeta))*(q_VT/q))*((S_v)/S)\
tau_r    = 
Cy_delr  = CL_aoav*dynpres_v*tau_r*br_frac*(S_v/s);
Cn_delr  = -CL_aoav*V_v*dynpres_v*tau_r*br_frac;

