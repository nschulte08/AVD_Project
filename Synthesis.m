clear; clc; close all;

%% Mission inputs
cruise_altitude = 51000; % ft
altitudes = [0, 30000, cruise_altitude]; % ft
M_cruise = 1.4;
num_pass = 19; % number of passengers
num_crew = 4; % number of crew members
range = 5468; % miles (supersonic cruise)
[~, ~, ~, ~, a] = ATMO(cruise_altitude, 'E');
V_cruise = M_cruise * a; % cruise velocity, ft/s 

%% Empirical inputs
SFC = 1.0;    % Empirical Placeholder for SFC based on Sadraey Table 4.6 for turbojet
LD_cruise = 9; % Cruise lift/drag from Fig 5.3 Nicolai

%% Interdisciplinary inputs
S = 15700; % wing area, ft^2
b = 120; % wingspan, m
lambda = 0.44; % wing taper ratio
M_max = 1.5;
CL_max = 2.6;

%% Phase-independent function calls
[weights, weight_fractions] = Weight_Buildup(num_pass, num_crew, convvel(V_cruise, 'ft/s', 'mph'), M_cruise, range, SFC, LD_cruise);
MTOW = weights.W_to; % Max takeoff weight, lbf
Vn_Diagram(MTOW, S, altitudes, M_cruise, M_max, CL_max);
[max_load, min_load] = Wing_Loading(b, convforce(MTOW, 'lbf', 'N'), lambda);

%% Takeoff
% Aerodynamics function calls
% Propulsion function calls
% Performance function calls
% Structures, W&B function calls
% Stability & control function calls

%% Climb
% Performance function calls

%% Cruise
% Aerodynamics function calls
% Propulsion function calls
% Performance function calls
% Stability & control function calls

%% Descent
% Performance function calls

%% Landing
% Aerodynamics function calls
% Propulsion function calls
% Performance function calls
% Structures, W&B function calls
% Stability & control function calls

