%{
    This function calculates the performance characteristics during the
    takeoff phase of the flight
---------------------------------------------------------------------------
Inputs:
ne:        Number of engines
V_TO:      Take off velocity (m/s)
T_TO:      Thrust at takeoff (N)
alt_to:    Airport runway altitude (m)
MTOW:      Aircraft maximum takeoff weight (N)
Sref:      Reference wing area (m^2)
CL_max:    Max lift at takeoff
CD0:       Zero lift drag coeff at takeoff
K:         drag K factor (1/piARe)
---------------------------------------------------------------------------
Outputs:
V_stall:   Stall velocity (m/s) 
V_TO:      Takeoff velocity (m/s)
S_G:       Ground roll distance (m)
S_TO:      Takeoff distance (total) (m)
BFL:       Balanced field length (m)
===========================================================================
%}
function [S_G, S_TO, BFL] = perf_takeoff(ne, V_stall, V_TO, T_TO, alt_to, MTOW, Sref, CL_TO, CD_TO)
% General Inputs:
g = 9.81;           % Gravity accel
OEI = 0;            % One engine out scenario, if yes [1] if no [0]
%--------------------------------------------------------------------------
if OEI == 0
    T_TO = T_TO;
elseif OEI == 1
    T_TO = ((ne-1)/ne)*T_TO; % [N] OEI scenario
end

%% ========================================================================
% Take off:
%--------------------------------------------------------------------------
[~, ~, ~, rho_TO, ~, ~, ~, ~, ~, SIGMA_TO] = ATMO(alt_to, 'M');

% ground roll:
mu = (0.02 + 0.3)/2;                 % friction coefficient [average value] (Yechout p.99)
q_TO = 0.5*rho_TO*V_TO^2;            % [N/m^2] dynamic pressure
%CL_TO = mu/(2*K);                   % optimum lift coefficient for T.O.
%CD_TO = CD0 + K*CL_opt^2;           % example 3.5 from Yechout (p.103)
D_TO = CD_TO*Sref*q_TO;              % [N] drag

if D_TO >= T_TO
	err_msg_1 = sprintf(' Drag at takeoff exceeds thrust at takeoff!! Cannot fly.... \n');
    err_msg_2 = sprintf('\n T_TO = %g [N] \n D_TO = %g [N]',T_TO, D_TO);
    err_msg = sprintf('%s %s',err_msg_1,err_msg_2);
	error('sytnax:requirements','\n\n %s \n\n',err_msg)
end

L_TO = CL_TO*Sref*q_TO;              % [N] lift
F_r = mu*(MTOW - L_TO);              % [N] rolling resistance
a_TO = (g/MTOW)*(T_TO - D_TO - F_r); % [m/s^2] average acceleration
S_G = V_TO^2/(2*a_TO);               % [m] ground distance
%--------------------------------------------------------------------------
% rotation:
t_R = 2.5;      % [s] rotation time
S_R = t_R*V_TO; % [m] distance covered during rotation
%--------------------------------------------------------------------------
% transition:
n = 1.2;                               % load factor
V_TR = 1.5*V_stall;                    % [m/s] transition velocity
Radius_TO = V_TR^2/(g*(n-1));          % [m] radius of transition
gamma_c = asind((T_TO - D_TO)/MTOW);   % [deg] climb angle
S_TR = Radius_TO*sind(gamma_c);        % [m] transition distance
h_TR = Radius_TO*(1 - cosd(gamma_c));  % [m] altitude gained during transition
%--------------------------------------------------------------------------
% climb:
h_obst = 35*0.3048; % [m] 35ft obstacle
if h_TR > h_obst
    S_C = 0; % [m] distance covered during climb
else
    S_C = (h_obst - h_TR)/tand(gamma_c);   % [m] distance covered during climb
end
%--------------------------------------------------------------------------
% total take-off distance:
S_TO = S_G + S_R + S_TR + S_C; % [m]
%--------------------------------------------------------------------------
% balanced field length: (Raymer eq. 17.110)
gamma_min = 0.027; % [deg]
G = (gamma_c - gamma_min)*pi/180;
U = 0.01*CL_TO + 0.02;
BFL = (0.863/(1+2.3*G))*((MTOW/Sref)/(rho_TO*g*CL_TO) + h_obst)*(1/((T_TO/MTOW) - U) + 2.7*0.3048) + (655*0.3048/sqrt(SIGMA_TO)); % [m] balanced field length
%--------------------------------------------------------------------------
fprintf('\n\n =========================== Take-off Results  =========================== \n');
fprintf('\n Max Take-off Gross Weight: MTOW = %g [N] = %g [lbf] \n', MTOW, MTOW*0.22480902);
fprintf('\n Stall Speed:    V_s  = %g [m/s] = %g [ft/s]', V_stall, V_stall*3.2808399);
fprintf('\n Take-off Speed: V_TO = %g [m/s] = %g [ft/s]', V_TO, V_TO*3.2808399);
fprintf('\n Obstacle clearance height: h_obst = %g [m] = %g [ft]', h_obst, h_obst*3.2808399);
fprintf('\n\n ------------------------------------------------------------- \n');
fprintf('\n Ground Roll:         S_G  = %g [m] = %g [ft]', S_G, S_G*3.2808399);
fprintf('\n Rotation Distance:   S_R  = %g [m] = %g [ft]', S_R, S_R*3.2808399);
fprintf('\n Transition Distance: S_TR = %g [m] = %g [ft]', S_TR, S_TR*3.2808399);
fprintf('\n Climb Distance:      S_C  = %g [m] = %g [ft]', S_C, S_C*3.2808399);
fprintf('\n\n Total Take-off Distance:      S_TO  = %g [m] = %g [ft]', S_TO, S_TO*3.2808399);
fprintf('\n\n Balanced Field Length:         BFL  = %g [m] = %g [ft]', BFL, BFL*3.2808399);
fprintf('\n\n ========================================================================= \n');

end
