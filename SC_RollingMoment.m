function [Cl_beta, Cl, ba_frac, ca_frac, Cl_dela, t2, phi1] = SC_RollingMoment(Lambda, TR, S, S_vt, S_vtpr, S_w, S_ref, b, c_r, d, sweep, AR, z_w, Cl_aoa, Cl_aoavt, rho)
%% Lateral Stability
%% Wing Geometry Defintions
% Lambda - Dihedral Angle
% TR - Taper Ratio
% S - Wing Area
% S_vt - Vertical Tail Area
% S_vtpr - V. Tail Area extended to fuselage
% S_w - Wing Area
% S_ref - Reference Area
% b - Wing Span
% c_r - Root Chord
% d - Max. Fuselage Depth
% sweep - Wing Sweep Angle
% AR - Aspect Ratio
% z_w - Distance b/w root chord and centerline along z-axis
%% Aerodynamic Coefficient Definitions
% Cl_aoa - Lift Coefficient due to angle of attack
% Cl_aoavt - Lift Coefficient due to angle of attack of the Vertical Tail

Cl_betabasic = 0; % Nicolai Figure 21.10 
Cl_betaD     = 0; % Nicolai Figure 21.10
z_v          = 0; % distance from mean aerodynamic chord of vertical stabilizer to Vertical CG Position
delL         = 0; % incremental change in the lift due to aileron deflection
I_xx         = 0; % Mass Moment of Inertia
y_D          = 0; % y-axis position of incremental drag (as averaged b/w surfaces(Wing and V.Tail))
phi2         = 0; % Required Bank Angle
C_Dr         = 0; % Drag Due to Roll
   
%% Common Terms
vwt         = (0.724+((3.06*(S_vtpr/S_ref))/(1+cos(sweep)))+(0.4*(z_w/d))+(0.009*AR)); % Nicolai 21.15

%% Stability Derivative Calculations
Cl_betaLam  = -0.25*Cl_aoa*Lambda*((2+(1+2*TR))/(3*(1+TR))); % Nicolai 21.12
Cl_betawing = Cl_betabasic+Cl_betaD+Cl_betaLam; % Nicolai 21.11
Cl_betavt   = (-Cl_aoavt*(S_vt/S_ref)*(z_v/b)*(vwt)); % Nicolai 21.14
Cl_beta     = Cl_betawing+Cl_betavt; % Nicolai 21.10
Cl_0        = 0; % Just a Placeholder

%% Roll Control Surface Sizing
ba_frac     = 0; % Aileron to wing span ratio
ca_frac     = 0; % Aileron to chord length ratio
tau_a       = ((1.5278*ca_frac^3)-(2.7083*ca_frac^2)+(2.2139*ca_frac)+0.0543);
y_i         = 0; % Aileron inboard y-position
y_o         = 0; % Aileron outboard y-position
in1         = (((y_o^2*(4*(TR-1)*y_o+3*b))-((4*(TR-1)*y_i^3)-(3*b*y_i^2)))/(6*b));
Cl_dela     = ((2*Cl_aoaw*tau_a*c_r)/(S*b))*(in1);
dela_max    = 0; % Maximum Aileron Deflection Angle [degrees]
Cl          = Cl_dela*dela_max;
L_A         = 2*delL*y_a;
P_ss        = sqrt((2*L_A)/(rho*(S_w+S_vt)*C_Dr*(y_D^3)));
phi1        = ((I_xx)/(rho*(y_D^3)*(S_w+S_vt)*C_Dr))*(log(P_ss^2));
P_dot       = (P_ss^2)/(2*phi1);
if phi1 > phi2
    t2      = sqrt((2*phi_des)/(P_dot));
    else 
    t2      = sqrt((2*phi1)/(P_dot))+((phi2-phi1)/(P_ss));
end
%% Plots
dela = -25:5:25;
beta = -4:2:24; 
for n = 1:length(dela)
    for j = 1:length(beta)
        BETA(j) = (beta(j)*(pi/180));
        Cl(n,j) = (Cl_0+(Cl_beta*BETA(j))+(Cl_dela*dela(n))+(Cl_dela*delr));
    end
end