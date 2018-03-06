%{
    Aerodynamic Function
    --------------------
    Accepts flight conditions and outputs aero coefficients

    Written by: Ben Holden
    Last Modified: 03/06/2018 (by Shawn McCullough)
===========================================================================
Inputs:
h = altitude [m]
M = Mach #
AR = aspect ratio
TR = taper ratio
Sref = referrence area [m^2]
SM = static margin
===========================================================================
Outputs:
M_crit = 
LD_max = 
CL = 
CD = 
CD0 = 
CM = 
CLa = 
CMa = 
sweep_rad = [radians]
alpha_trim = [radians]
===========================================================================
%}
function [M_crit, LD_max, CL, CD, CD0, CM, CLa, CMa, sweep_rad, alpha_trim] = aerofunky3(h, M, AR, TR, Sref, SM)
%% other known inputs
tcmax = 0.16; % known because of ... (this is internal to aero team... know one else uses it)
Swet = 2.1*Sref; % [m^2] Placeholder, Kinda a guess, 2*Sref for flying wing (no one else uses this... it is internal to aero team)
e = 0.85;           % Oswald
K = 1/(pi*AR*e);	% Drag K factor

%cbar = 1.0;        % [m] Placeholder for mean geometric chord
c_r = 1;            % [m] Placeholder for root chord
cbar = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));  % Nicolai Pg 580

% Calculate sweep:
M_perp = 0.7; % perpendicular Mach #
sweep_deg = acosd(M_perp/M);
sweep_rad = sweep_deg*pi/180; % Convert to radians for formulas

% sweep = f(M_cr); Mcrit = 
% We also need to calculate b_eff and AR_eff with the sweep angle

%% Lift Calculations

%--------Ideal Wing Lift Coefficient, Saedray Ch5--------%
% W_cs = 1.0;                                %Placeholder %Weight at start of cruise
% W_cf = 1.0;                                %Placeholder %Weight at end of cruise
%[~, a, ~, rho] = atmosisa(h);      % Cruise Atmospheric Properties
[~, ~, ~, rho, a, mu, ~, ~, ~, ~] = ATMO(h, 'M', 'N');
V = M*a;                           % Cruise Speed (m/s)
%Calculate CL for steady level flight
%CL_cruise = 2*Weight/(rho*V^2*Sref);     %Ideal wing lift coefficient for cruise flight

%---------- Wing Lift Curve Slopes ----------%
%Subsonic 
Cla = 1.8*pi*(1 + 0.8*tcmax); %Airfoil Cl_alpha, Saedray pg 179  
Beta_sub = sqrt(abs(1-M^2));
nu = Cla/(2*pi/Beta_sub);
CLa_sub = 2*pi*AR/(2+sqrt(4+((AR^2*Beta_sub^2)/nu^2)*(1+tan(sweep_rad)^2/Beta_sub^2))); 

%Supersonic 
%CLa_super = 4/(sqrt(M_cr^2-1));
CLa_super = 4/(sqrt(M^2-1))*(1-1/(2*AR*sqrt(M^2-1))); %Straight,tapered wings in supersonic flow

% Calculate Angle of Attack for Trim
CM0 = 1.0;                  %Placeholder, CM0 comes from airfoil
CMa_sub = -CLa_sub*SM;
CMa_super = -CLa_super*SM;

alpha_trim_sub = -CM0/CMa_sub;     % radians
alpha_trim_super = -CM0/CMa_super; % radians

%Calculate lift coefficients based on CLa and alpha
CL_sub = alpha_trim_sub*CLa_sub;                     % Alpha needs to be in radians here
CL_super = alpha_trim_super*CLa_super;               % Alpha needs to be in radians here
CM_sub = CMa_sub*alpha_trim_sub;
CM_super = CMa_super*alpha_trim_super;

%% Transonic Wing Aerodynamics

Mach_dd = 0.95/cos(sweep_rad) - (tcmax)/cos(sweep_rad)^2 - CL_sub/(10*cos(sweep_rad)^3);        %Drag divergence Mach number
M_crit = Mach_dd - (0.1/80)^(1/3);                                                              %Critical Mach number

%% Drag Analysis

% Subsonic Drag
% These are inputs from Roskam ch4 Figures:
%mew = 1.422*10^-5;                             % Dyn Viscocity for 15,000 m
Re = rho*V*cbar/mu;                             % Reynolds Number
Cfw = 0.455/(log10(Re)^2.58);                   % Turbulent flat plate friction coefficient of wing

%Alternate method:Subsonic CD0 is estimated from skin friction
CD0_sub = 1.2*Cfw*Swet/Sref;                            

%Calculate subsonic, steady level LD max
LD_max = 0.5*sqrt((pi*AR*e)/CD0_sub);                 

% Calculate Lift-induced drag as a f(twist,CL)
% Zero unless twist is nonzero, then from Roskam ch 4 find parameters
eta_t = 0;           
v_t = 0;
w_t = 0;
CDi_sub = CL_sub^2*K + 2*pi*CL_sub*eta_t*v_t + 4*pi^2*eta_t^2*w_t;  

% Total Subsonic CD
CD_sub = CD0_sub + CDi_sub;    

% Supersonic Drag
Cfw_super = Cfw/(1 + .144*M^2)^0.65;                                %Corrected Skin Friction Coefficient for Compressibility
CDf_super = Cfw_super*Swet/Sref;                                    %Supersonic friction drag coefficient

clms = 1.0;                                                         %Placeholder %Camber Line Mean Square
tdms = 1.0;                                                         %Placeholder %Thickness distribution mean square
CD_wave = 4*alpha_trim_super/(M^2 - 1)^0.5*(alpha_trim_super^2 + clms + tdms); %Section Wave drag coefficient

CD0_super = CDf_super + CD_wave;                                    %Supersonic zero-lift drag coefficient

Kw = .175;                                                          %Wing wave drag factor, Howe eq 6.18
CDi_super = (0.24/AR + Kw*(M^2-1)^.5)*CL_super^2;                   %Supersonic induced drag coefficient (Howe)

% Total Supersonic CD
CD_super = CD0_super + CDi_super;                                   %Total supersonic CD:Sum of zero-lift + lift-induced drag

%% Aero Coefficient Assignments
if M < 0.8
    CL = CL_sub;
    CD = CD_sub;
    CD0 = CD0_sub;
    CLa = CLa_sub;
    alpha_trim = alpha_trim_sub;
    CM = CM_sub;
    CMa = CMa_sub;
% -------------------------------------------------------------------------
elseif M > 1.1
    CL = CL_super;
    CD = CD_super;
    CD0 = CD0_super;
    CLa = CLa_super;
    alpha_trim = alpha_trim_super;
    CM = CM_super;
    CMa = CMa_super;
% -------------------------------------------------------------------------
else % (transonic)
    err_msg = sprintf(' Cannot calculate transic aerodynamics. \n\n   M < 0.8  OR  M > 1.1 ');
	error('sytnax:requirements','\n\n %s \n\n',err_msg)
    % using linear interpolation between to end points: (doesn't work --> complex numbers...?)
    %{
    CL_data    = [CL_sub, CL_super];
    CD_data    = [CD_sub, CD_super];
    CD0_data   = [CD0_sub, CD0_super];
    CM_data    = [CM_sub, CM_super];
	CLa_data   = [CLa_sub, CLa_super];
    CMa_data   = [CMa_sub, CMa_super];
    Mach       = [0.8, 1.1];
    
    CL  = interp1(Mach, CL_data, M);
    CD  = interp1(Mach, CD_data, M);
    CD0 = interp1(Mach, CD0_data, M);
    CLa = interp1(Mach, CM_data, M);
    CM  = interp1(Mach, CLa_data, M);
    CMa = interp1(Mach, CMa_data, M);
    alpha_trim = nan;
    %}
end

end