%This script calculates the aerodynamic coefficients based on flight
%conditions

%Required inputs: Cruise condition (Mach, altitude, Weight Fractions), t/cmax, 
%Sref, AR, Swet, sweep angle

Mc = 1.0;                               %Placeholder, cruise Mach number
AR = 1.0;                               %Placeholder, AR needs to be f(sweep)
K = 1/(pi*AR*e);                        %Drag K factor
tcmax = 1.0;                            %Placeholder %Max t/c ratio 
altcr = 1.0;                            %Cruise Altitude (units?)
Sref = 1.0;                             %Wing Area, (units?)

%% Airfoil Cl_alpha

Cla = 1.8*pi*(1 + 0.8*tcmax);         %Airfoil Cl_a, Saedray pg 179

%% Wing Lift Coefficient Calculation - Saedray Ch 5

Wi = 1.0;                               %Placeholder %Weight at start of cruise
Wf = 1.0;                               %Placeholder %Weight at end of cruise
Wavg = 0.5*(Wi + Wf);                   %Average Cruise Weight
[Tcr, acr, Pcr, rhocr] = ATMO(altcr);   %Cruise Atmospheric Properties
Vc = 1.0;                               %Cruise Speed (units?)
CL_ideal = 2*Wavg/(rhocr*Vc^2*Sref);    %Ideal wing lift coefficient for cruise flight

%% Subsonic Wing Aerodynamics

%Subsonic CL_alpha
Bsub = sqrt(1-Mc^2);
nu = Cla/(2*pi/Bsub);
CLa_sub = 2*pi*AR/(2+sqrt(4+((AR^2*Bsub^2)/nu^2)*(1+tan(sweep)^2/Bsub^2)));            %Subsonic wing lift curve slope

%Subsonic Drag
%These are inputs from Roskam ch4 Figures:
Rew = rhocr*Vc*cbar/mu;                             %Reynolds Number for wing at cruise
Cfw = 0.455/(log10(Rew)^2.58);                      %Turbulent flat plate friction coefficient of wing

CD0 = 1.2*Cfw*Swet/Sref;                            %Alternate method:Subsonic CD0 is estimated from skin friction
CLCD_max = 0.5*sqrt((pi*AR*e)/CD0);                 %Subsonic steady level CL/CD Max

%Zero unless twist is nonzero, then from Roskam ch 4 find parameters
eta_t = 0;           
v_t = 0;
w_t = 0;

CDL = CL^2*K + 2*pi*CL*eta_t*v_t + 4*pi^2*eta_t^2*w_t;  %Lift-induced drag = f(twist,CL)

%Total Subsonic CD
CD = CD0 + CDL;                                     %Subsonic drag coefficient


%% Transonic Wing Aerodynamics

Mach_dd = 0.95/cos(sweep) - (t/c)/cos(sweep)^2 - CL/(10*cos(sweep)^3);      %Drag divergence Mach number
M_cr = Mach_dd - (0.1/80)^(1/3);                                            %Critical Mach number


%% Supersonic Wing Aerodynamics

%Supersonic CLa
CLa_super = 4/(sqrt(Mc^2-1));

%Supersonic Drag
Cfw_super = Cfw/(1 + .144*Mc^2)^0.65;                       %Corrected Skin Friction Coefficient for Compressibility
CDf_super = Cfw_super*Swet/Sref;                            %Supersonic friction drag coefficient

clms = 1.0;                         %Placeholder %Camber Line Mean Square
tdms = 1.0;                         %Placeholder %Thickness distribution mean square
Cd_wave = 4*alpha/(Mc^2 - 1)^0.5*(alpha^2 + clms + tdms;    %Section Wave drag coefficient

CD0_super = CDf_super + CD_wave;                            %Supersonic zero-lift drag coefficient

Kw = .175;                                                  %Wing wave drag factor, Howe eq 6.18
CDL_super = (0.24/AR + Kw*(Mc^2-1)^.5)*CL^2;                %Supersonic induced drag coefficient (Howe)

%Total Supersonic CD
CD_super = CD0_super + CDL_super;                           %Total supersonic CD:Sum of zero-lift + lift-induced drag


