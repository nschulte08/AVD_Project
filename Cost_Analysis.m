% This script is built for the cost analysis
clear; close all; clc;
%--------------------------------------------------------------------------
W_e = 49800.0;          %Empty Weight [in lbs]
V_max = 868.976;        %Max Speed [eq's need in knots]
Qp = 300.0;             %Quantity produced
Qd = 8;                 %Placeholder, quanity produced for RTDE phase
Q = Qp + Qd;
T_sl = 51000.0;         %Comes from propulsion
M_max = 1.5;            %Placeholder, comes from performance
Temp_Inlet = 3200.0;	%Placeholder, comes from propulsion
CPI_equiv = 29/20;      %Equivalence of today's dollar value
%--------------------------------------------------------------------------
%Airframe Engineering Cost-Part of DT&E
AE_hr = 4.86*(W_e^0.777)*(V_max^0.894)*(Q^0.163);	%E is engineering hours, We is empty wt in lbs, Vmax in knots, Qp is quantity produced
AE_rate=2.576*(2017)-5058;                          %Eq from Fig 24.4 Nicolai
AE_cost = AE_hr*AE_rate;                            %Cost of airframe engineering = E*Engineering Rate
%--------------------------------------------------------------------------
%Development Support Cost (D)
DSC = (66*(W_e^0.63)*(V_max^1.3))*CPI_equiv;        %Cost for development support, We in lbs, Vmax in knots
%--------------------------------------------------------------------------
%Flight Test Op Cost (F)
FTC = (1852*(W_e^0.325)*(V_max^0.822)*(Qd^1.21))*CPI_equiv;	%Cost for flight test, same units as above
%--------------------------------------------------------------------------
%Tooling costs
T_hrd = 5.99*(W_e^.777)*(V_max^.696)*(Q^.263);	%Tooling hours, same units as above
T_rate = 2.883*(2017) - 5666;   
T_cost = T_hrd*T_rate;
%--------------------------------------------------------------------------
%Manufacturing Cost
L_hr = 7.37*(W_e^.82)*(V_max^.484)*(Q^.641);	%Labor hours, same units as above
L_rate = 2.316*(2017)-4552;
L_cost = L_hr*L_rate;
%--------------------------------------------------------------------------
%Quality Control Cost
QC_hr = 0.13*L_hr;	%Quality control hours
QC_rate = 2.6*(2017) -5112;
QC_cost = QC_hr*QC_rate;
%--------------------------------------------------------------------------
%Material Costs
M_cost = (16.39*(W_e^.921)*(V_max^.621)*(Q^.799))*CPI_equiv;
%--------------------------------------------------------------------------
%Propulsion & Avionics Cost
P_per = 3*2306*(0.043*T_sl + 243.3*M_max +.969*Temp_Inlet - 2228)*CPI_equiv; %TSL is max thrust in lbs, Temp is in Rankine
P_cost = Qd*P_per;
%--------------------------------------------------------------------------
RTDE_Ac_Cost = AE_cost + DSC + FTC + T_cost + L_cost + QC_cost + M_cost + P_cost;
%--------------------------------------------------------------------------
disp(num2str(RTDE_Ac_Cost,'RTDE = %.2f'))
disp(num2str(AE_cost,'Airframe Engineering = %.2f'))
disp(num2str(DSC,'Development Support = %.2f'))
disp(num2str(FTC,'Flight Test = %.2f'))
disp(num2str(T_cost,'Tooling = %.2f'))
disp(num2str(L_cost,'Manufacturing = %.2f'))
disp(num2str(QC_cost,'Quality Control = %.2f'))
disp(num2str(M_cost,'Material = %.2f'))
disp(num2str(P_cost,'Propulsion/Avionics = %.2f'))