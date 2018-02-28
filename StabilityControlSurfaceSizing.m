%% Control Surface Sizing==================================================
    %By: Christian Allen
    %12/19/18
%% Flight Phases===========================================================
%Class()
%Phase(B)
%% Aircraft Geometry/Variables=============================================
S        
TR       
c        
c_root   
b       
b_vt   
rho    
CL_aoaw
delL%Incremental change in the lift due to airleron deflection   
y_a%CG to aileron center location y-position 
S_w%Wing Area    
S_vt%Vertical Tail Area
C_Dr%Drag Due to Roll
y_D%Incremental drag y_position from CG point    
I_xx%Mass moment of inertia about x-axis    
phi2     
phi_req  
Cmac_wf 
mac
T 
D 
miu%Ground Friction Coefficient 
N%Normal Force 
m%mass 
W%Weight 
x_mg%Main gear x-position
x_cg%Center of Gravity x-position
z_D%Incremental drag z_position from CG point 
z_mg%Main gear z-position 
z_T%Thrust vector z-position 
L_wf%Lift due to wing-fuesl 
x_acwf%Aerodynamic Center x-position
a%accelteration 
z_cg%Center of Gravity z-position 
I_yymg%Mass moment of inertia about the main gear 
theta_ddot
x_ach%Aerodynamic Center x-position of H.tail 
S_h%Horizontal Tail Area
CL_aoah
aoa_h
dynpres_h%H. Tail dynamic pressure ratio
V_h
q_bar
mac
CL_aoa
CL_1%Steady-State aircraft lift coefficient at cruising flight
CL_0
Cm_aoa
T_L%operative Engine Thrust  
y_T%Engine y_position from centerline   
CL_aoav =2*pi
dynpres_v%V. Tail dynamic pressure ratio
b_r      
U1
CD_y%Side Drag Coefficient 
K_f1 
l_vt%Vertical Tail Arm
S_v%Vertical Tail Area  
K_f2 
%% Aileron=================================================================
ba_frac   = 0.2
ca_frac   = 0.25 
tau_a     = 0.48
y_i       = 0.8
y_o       = 1
in1       = (((y_o^2*(4*(TR-1)*y_o+3*b))-((4*(TR-1)*y_i^3)-(3*b*y_i^2)))/(6*b));
Cl_dela   = ((2*CL_aoaw*tau_a*c_root)/(S*b))*(in1);
dela_max = 25;
Cl        = Cl_dela*dela_max;
L_A       = 2*delL*y_a;
P_ss      = sqrt((2*L_A)/(rho*(S_w+S_vt)*C_Dr*(y_D^3)));
phi1      = ((I_xx)/(rho*(y_D^3)*(S_w+S_vt)*C_Dr))*(log(P_ss^2));
P_dot     = (P_ss^2)/(2*phi1);
if phi1 > phi_req
t2        = sqrt((2*phi_des)/(P_dot));
    else  phi1 < phi_req
t2        = sqrt((2*phi1)/(P_dot))+((phi2-phi1)/(P_ss));
end
%% Elevators(Takeoff)======================================================
be_frac   = 1
dele_max  = 30
L_wf      = 0.5*rho*(V_r^2)*CL_TO*S_ref;
D_TO      = 0.5*rho*(V_r^2)*CD_TO*S_ref;
Mac_wf    = 0.5*rho*(V_r^2)*Cmac_wf*S_ref*mac;
a_TO      = (T-D-(miu*N))/m;
M_w       = (W*(x_mg-x_cg));
M_D       = (D*(z_D-z_mg));
M_T       = (T*(z_T-z_mg));
ML_wf     = (L_wf*(x_mg-x_acwf));
M_a       = (m*a*(z_cg-z_mg));
L_h       = ((ML_wf+Mac_wf+M_a-M_w+M_D-M_T-(I_yymg*theta_ddot))/(x_ach-x_mg))
CL_h      = ((2*L_h)/(rho*(V_r^2)*S_h));
tau_e     = (((CL_h/CL_aoah)-aoa_h)/dele_max);
ce_frac   = 0.5
Cm_dele   = -CL_aoah*dynpres_h*V_h*be_frac*tau_e;
CL_dele   = CL_aoah*dynpres_h*(S_h/S)*be_frac*tau_e;
CL_hdele  = CL_aoah*tau_e;
dele_trim = -((((((T*z_T)/(q_bar*S*mac))*CL_aoa)+(CL_1-CL_0)*Cm_aoa))/((CL_aoa*Cm_dele)-(Cm_aoa*CL_dele)))
%% Rudder==================================================================
%Asymmetric Thrust=========================================================
br_frac   = 1
delr_max  = 30;
V_conmin  = 16 
N_A       = -T_L*y_T;
Cn_delr   = N_A/(q*S*b*delr_max);
tau_r     = (Cn_delr*b_vt/(-CL_aoav*V_v*dynpres_v*b_r));
cr_frac   = 0.5
%Cross-wind Landing========================================================
V_v       = ((l_vt*S_v)/(b*S));
V_w       = 28
V_t       = sqrt((U1^2)+(V_w^2));
S_s       = 84
F_w       = 0.5*rho*(V_w^2)*S_s*CD_y;
b_frac    = 1
c_frac    = 1
beta      = atan(V_w/U1);
Cn_beta   = K_f1*CL_aoav*(vwt)*((l_vt*S_v)/S)
Cy_beta   = K_f2*CL_aoav*(vwt)*((S_v)/S)
tau_r     = 0.5
Cy_delr   = CL_aoav*dynpres_v*tau_r*br_frac*(S_v/S);
Cn_delr   = -CL_aoav*V_v*dynpres_v*tau_r*br_frac; 

