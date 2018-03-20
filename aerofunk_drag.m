function [CD, CD0] = aerofunk_drag(h, M, Sref, CL, alpha, AR, TR)
%This function returns the drag coefficients based on CL and flight phase

%Inputs:
%h:     Altitude, m
%M:     Mach number for flight phase
%Sref:  Wing Reference Area, m^2
%CL:    Lift coefficient based on flight profile
%AR:    Aspect Ratio (it is being calculated based on sweep in synthesis
%TR:    Taper Ratio


%Outputs:
%CD:    Drag Coefficient (will output based on input Mach number)
%CD0:   Zero-lift drag coefficient (will output based on input Mach number)

% ==================================================================================== %

% Preliminary Calculations:
tmax = 2.3;         % Maximum Thickness
tcmax = 0.16;       % T/c max; variable to iterate (this is internal to aero team... ow one else uses it)
c_r = tmax/tcmax;            % [m] Root chord = max thickness / tcmax ratio
cbar = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));  % Nicolai Pg 580

Swet = 2.1*Sref;    % [m^2] Placeholder, Kinda a guess, 2*Sref for flying wing (no one else uses this... it is internal to aero team)
e = 0.85;           % Oswald
K = 1/(pi*AR*e);	% Drag K factor

% Calculate sweep:
% This has to be limited to 70 ish degrees!
M_perp = 0.7; % perpendicular Mach #, variable to iterate
sweep_deg = acosd(M_perp/M);

if sweep_deg > 70                       %Limit the sweep angle
    sweep_deg = 70;
end

sweep_rad = sweep_deg*pi/180; % Convert to radians for formulas

% Get atmospheric properties
[~, ~, ~, rho, son, mu, ~, ~, ~, ~] = ATMO(h, 'M');
V = son * M;            %Velocity based on Mach number input
Re = rho*V*cbar/mu;                             % Reynolds Number
Cfw = 0.455/(log10(Re)^2.58);                   % Turbulent flat plate friction coefficient of wing


if M < .9
    %---------------------------------- Subsonic Drag ----------------------------------%

    % Subsonic CD0 is estimated from skin friction
    CD0_sub = 1.2*Cfw*Swet/Sref;                            

    % Calculate Lift-induced drag as a f(twist,CL)
    % Zero unless twist is nonzero, then from Roskam ch 4 find parameters
    eta_t = 0;           
    v_t = 0;
    w_t = 0;
    CDi_sub = CL^2*K + 2*pi*CL*eta_t*v_t + 4*pi^2*eta_t^2*w_t;  

    % Total Subsonic CD
    CD0 = CD0_sub;
    CD = CD0 + CDi_sub;    

    
elseif M > 1.1
    % ------------------------- Supersonic Drag ------------------------- %
    Cfw_super = Cfw/(1 + .144 * M^2)^0.65;                                %Corrected Skin Friction Coefficient for Compressibility
    CDf_super = Cfw_super * Swet/Sref;                                    %Supersonic friction drag coefficient

    clms = 1.0;                                                         %Placeholder %Camber Line Mean Square (pg 65 Nicolai)
    tdms = 1.0;                                                         %Placeholder %Thickness distribution mean square (pg 65 Nicolai)
    CD_wave = 4 * alpha/(M^2 - 1)^0.5 * (alpha^2 + clms + tdms); %Section Wave drag coefficient

    CD0_super = CDf_super + CD_wave;                                    %Supersonic zero-lift drag coefficient

    Kw = .175;                                                          %Wing wave drag factor, Howe eq 6.18
    CDi_super = CL^2 * (0.24/AR + Kw * (M^2-1)^.5);                   %Supersonic induced drag coefficient (Howe)

    % Total Supersonic CD
    CD0 = CD0_super;
    CD = CD0 + CDi_super;                                   %Total supersonic CD:Sum of zero-lift + lift-induced drag
    
end
