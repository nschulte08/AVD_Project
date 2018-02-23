function [ M_crit, CL_sub, CD_sub, LD_max, CD0_sub, CL_super, CD_super, CD0_super, CMa, sweep] = aerofunky(alpha, M_cr, AR, alt_cr, Sref, W_cs, W_cf )

% Accepts flight conditions and outputs aero coefficients
% Required inputs: Cruise condition (Mach, altitude, Weight Fractions), t/cmax, 
% Sref, AR, Swet (maybe)

Swet = 1.05*Sref;                       %Placeholder, Kinda a guess
e = 0.85;
K = 1/(pi*AR*e);                        %Drag K factor
tcmax = 1.0;                            %Placeholder %Max t/c ratio 
SM = .15;                               %Static Margin as % of mean geom chord
cbar = 1.0;                             %Placeholder for mean geometric chord

%Calculate sweep: sweep = f(M_cr); Mcrit = 
sweep = 1.0;                            %Placeholder

%We also need to calculate b_eff and AR_eff with the sweep angle

%% Wing Lift Coefficient Calculation - Saedray Ch 5

%W_cs = 1.0;                                %Placeholder %Weight at start of cruise
%W_cf = 1.0;                                %Placeholder %Weight at end of cruise
Wavg = 0.5*(W_cs + W_cf);                   %Average Cruise Weight
[~, acr, ~, rhocr] = atmosisa(alt_cr);      %Cruise Atmospheric Properties
V_cr = M_cr*acr;                            %Cruise Speed (units?)

%Calculate CL for steady, level cruise
CL_cruise = 2*Wavg/(rhocr*V_cr^2*Sref);     %Ideal wing lift coefficient for cruise flight

%% Subsonic Wing Aerodynamics

%Subsonic wing lift curve slope
Cla = 1.8*pi*(1 + 0.8*tcmax);                       %Airfoil Cl_alpha, Saedray pg 179  
Beta_sub = sqrt(1-M_cr^2);
nu = Cla/(2*pi/Beta_sub);

CLa_sub = 2*pi*AR/(2+sqrt(4+((AR^2*Beta_sub^2)/nu^2)*(1+tan(sweep)^2/Beta_sub^2))); 

%Subsonic lift coefficient
CL_sub = alpha*CLa_sub;

%Subsonic Drag
%These are inputs from Roskam ch4 Figures:
mew = 1.0;                                          %Placeholder for mu
Rew = rhocr*V_cr*cbar/mew;                          %Reynolds Number for wing at cruise
Cfw = 0.455/(log10(Rew)^2.58);                      %Turbulent flat plate friction coefficient of wing

%Alternate method:Subsonic CD0 is estimated from skin friction
CD0_sub = 1.2*Cfw*Swet/Sref;                            

%Calculate subsonic, steady level LD max
LD_max = 0.5*sqrt((pi*AR*e)/CD0_sub);                 

%Zero unless twist is nonzero, then from Roskam ch 4 find parameters
eta_t = 0;           
v_t = 0;
w_t = 0;

%Calculate Lift-induced drag as a f(twist,CL)
CDi_sub = CL_sub^2*K + 2*pi*CL_sub*eta_t*v_t + 4*pi^2*eta_t^2*w_t;  

%Total Subsonic CD
CD_sub = CD0_sub + CDi_sub;                                     

%% Transonic Wing Aerodynamics

Mach_dd = 0.95/cos(sweep) - (tcmax)/cos(sweep)^2 - CL_sub/(10*cos(sweep)^3);        %Drag divergence Mach number
M_crit = Mach_dd - (0.1/80)^(1/3);                                                  %Critical Mach number

%% Supersonic Wing Aerodynamics

%Supersonic CLa
CLa_super = 4/(sqrt(M_cr^2-1));
%CLa_super = 4/(sqrt(M_cr^2-1))*(1-1/(2*AR*sqrt(M_cr^2-1));     %Straight,tapered wings in supersonic flow
CL_super = alpha*CLa_super;

%Supersonic Drag
Cfw_super = Cfw/(1 + .144*M_cr^2)^0.65;                       %Corrected Skin Friction Coefficient for Compressibility
CDf_super = Cfw_super*Swet/Sref;                            %Supersonic friction drag coefficient

clms = 1.0;                         %Placeholder %Camber Line Mean Square
tdms = 1.0;                         %Placeholder %Thickness distribution mean square
CD_wave = 4*alpha/(M_cr^2 - 1)^0.5*(alpha^2 + clms + tdms);    %Section Wave drag coefficient

CD0_super = CDf_super + CD_wave;                            %Supersonic zero-lift drag coefficient

Kw = .175;                                                  %Wing wave drag factor, Howe eq 6.18
CDi_super = (0.24/AR + Kw*(M_cr^2-1)^.5)*CL_super^2;                %Supersonic induced drag coefficient (Howe)

%Total Supersonic CD
CD_super = CD0_super + CDi_super;                           %Total supersonic CD:Sum of zero-lift + lift-induced drag


%% Pitching Moment Calculation
if M_cr > 1
    CLa = CLa_super;
else
    CLa = CLa_sub;
end

CMa = -CLa*SM;

%Do we want to have if M_cr < 1
            %CD = CD_sub and
            %if M_cr > 1 CD=CD_super etc or just return them all?

end

