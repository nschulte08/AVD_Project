function [ DOC_per_km ] = DOC_SSBJ(MTOW, W_empty, W_fuel, V_cr, Tto, ne, R, TOF, PAX)
% This script calculates the operating cost for the SSBJ 

%From synthesis script: 
% MTOW       Input in N
% W_empty	Input in N
% W_fuel     Input in N
% V_cr       Input as m/s
% Tto        Total thrust in kN
% ne         Number of engines
% R          Total range (km)
% TOF        TOF in hr
% PAX        PAX / number of seats

% Inputs:
c = 4;                  % Number of preproduction aircraft
Dac = 8.255;            % Aircraft deflator
De = 8.02;              % Engine deflator
Dl = 6.624;             % Labor deflator
Df = 25;                 % Fuel price deflator
nb = 150;               % PlaceholderBreakeven number (150?)
sew = .15;              % Given as between .1 and .25

%% Research & Development Costs:
V_c = convvel(V_cr, 'm/s','mph');

V = V_c/1000;           % Cruise speed / 1000   (Mm/hr)
maf = W_empty/9.81/1000;       % Airframe mass in 1000*kg
CDaf = Dac * 6.2e6 * maf;       % Airframe development cost Maybe multiply by .82V^2 ? for supersonic

Caf = nb^-0.322 * CDaf / c;     % Cost per aircraft
Ce = De * 24000 * exp(-7*sew) * Tto;    % Per engine cost

Ct = Caf + ne * Ce;         % Total aircraft cost

%% Operating Costs:
%Inputs:
% For Cakm1
tloss = 0.3;            % Average time loss (hr)
tflight = R/(1000 * V * 0.8) + 0.3 * V; %Use V calculated in RD cost part, not V_cr !
Vb = R / (tloss + tflight);     % Block speed (km/h)
mto = MTOW/9.81/1000;           % Takeoff mass in 1000*kg

% For Cakm2
mftrip = (W_fuel/1.05)/9.81/1000;   % Fuel used during flight exluding reserve (1000*kg)
tb = 8;                         % Placeholder Block time (from TOF?)
U = 2650 + 2100 * (1 - exp(-tb/2 + 0.5));   % Utilization

% For Cakm3
IR = .01;             % Placeholder, Insurance Rate (prosposed as 1%)

% For Cakm4
Da = 15;            % Placeholder, Depreciation period

% For Cakm5
AFCL = (2.14 + .0079*mto + .0046*PAX)/sqrt(V);      % Number of labor hours per airframe per flight cycle
AFHL = (3.08 + 0.032*maf + 0.0041*PAX)/sqrt(V);     % Number of labor hours per airframe per flight hour

% For Cakm6
TET =  2000;                                        % Turbine Entry Temp (K) per Rustin
mtbr = 3604 * tb^0.28 * exp(-0.000324*TET);          % Mean time between repairs (hr)
meng = sew*Tto/9.81;                                % Mass of engine in 1000*kg
EL = 0.143 + (1452 + 530 * meng)/mtbr;              % # of labor hours per engine per flight hour (hr)

% For Cakm7
AFCM = (7.23 + 0.096*mto + 0.020*PAX)* .82*V^2;      % Materials costs per airframe per flight cycle
AFHM = (6.51 + 0.028*maf + 0.025*PAX)* .82*V^2;      % Materials costs per airframe per flight hour

% DOC Calculations
Cakm1 = Dl * (5.66 * mto^0.7) / Vb;                     % Cost of flight crew ($/km)
Cakm2 = Df * 1.02 * (mftrip * 40 + ne * tb * .145)/R;   % Fuel cost, ($/km), G650 is 2840 $/hr
Cakm3 = IR * (Caf + ne*Ce)/(U * Vb);                    % Cost of Insurance ($/km)
Cakm4 = 0.9 * (1.04*Caf + ne*Ce*1.3) / (Da*U*Vb);       % Cost of Depreciation ($/km)
Cakm5 = Dl * 4.82 * (AFCL/R + AFHL/Vb);                 % Airframe labor costs
Cakm6 = (Dl * 4.82 * ne * EL/Vb)*.82*V^2;               % Engine labor costs
Cakm7 = Dac * (AFHM/Vb + AFCM/R);                       % Airframe parts costs
Cakm8 = ne * (De * 0.626 + 0.045 * Ce/mtbr)/Vb;         % Engine parts costs
Cakm9 = 1.4 * (Cakm5 + Cakm6);                          % Maintenance burden, indirect maintenance costs

%loadfactor = 0.8;                      % Average for most airlines

% Calculate total cost per km:
DOC_per_km = (Cakm1 + Cakm2 + Cakm3 + Cakm4 + Cakm5 + Cakm6 ...
    + Cakm7 + Cakm8 + Cakm9);

% Calculate total cost per km per PAX
%TDOC = DOC_per_km/(loadfactor*PAX);

% Just for comparison (used in validation)
G650_per_km = 29977/convlength(4000, 'mi', 'km');

end

