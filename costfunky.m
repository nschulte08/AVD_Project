%{
    Calculate the RDTE and DOC costs
---------------------------------------------------------------------------
INPUTS:
MTOW        = max take off weight [N]
W_empty     = empty weight [N]
W_fuel      = fuel weight [N]
V_cr        = cruise speed [N]
Tto         = Total thrust [kN] ????????????????????????????????
ne          = Number of engines
Range       = Total range [km]
TOF         = TOF [hr]
PAX         = number of passengers
T_SL_max    = max sea levell thrust [N]
---------------------------------------------------------------------------
OUTPUTS:
RTDE_Cost = research developement and E cost [$]
DOC_Cost  =  
---------------------------------------------------------------------------
%}
function [ RTDE_Cost, DOC_Cost ] = costfunky(MTOW, W_empty, W_fuel, V_cr, ne, Range, PAX, T_SL_max, M_max)
%--------------------------------------------------------------------------
% convert units for RDTE cost analysis:
W_e = convforce(W_empty,'N','lbf');
%W_f = convforce(W_fuel,'N','lbf');
T_sl = convforce(T_SL_max,'N','lbf');
V_max = convvel(V_cr, 'm/s','kts');
%--------------------------------------------------------------------------


%% ========================================================================
% RDTE:
%--------------------------------------------------------------------------
Qp = 300.0;             % Quantity produced
Qd = 8;                 % Placeholder, quanity produced for RTDE phase
Q = Qp + Qd;
Temp_Inlet = 3200.0;	% Placeholder, comes from propulsion
CPI_equiv = 29/20;      % Equivalence of today's dollar value
%--------------------------------------------------------------------------
% Airframe Engineering Cost-Part of DT&E
AE_hr = 4.86*(W_e^0.777)*(V_max^0.894)*(Q^0.163);	% E is engineering hours, We is empty wt in lbs, Vmax in knots, Qp is quantity produced
AE_rate=2.576*(2018)-5058;                          % Eq from Fig 24.4 Nicolai
AE_cost = AE_hr*AE_rate;                            % Cost of airframe engineering = E*Engineering Rate
%fprintf('\n Cost of airframe engineering = $%g ', AE_cost);
%--------------------------------------------------------------------------
% Development Support Cost (D)
DSC = (66*(W_e^0.63)*(V_max^1.3))*CPI_equiv;        %Cost for development support, We in lbs, Vmax in knots
%fprintf('\n Cost for development support = $%g ', DSC);
%--------------------------------------------------------------------------
% Flight Test Op Cost (F)
FTC = (1852*(W_e^0.325)*(V_max^0.822)*(Qd^1.21))*CPI_equiv;	%Cost for flight test, same units as above
%fprintf('\n Cost for flight test = $%g ', FTC);
%--------------------------------------------------------------------------
%Tooling costs
T_hrd = 5.99*(W_e^.777)*(V_max^.696)*(Q^.263);	%Tooling hours, same units as above
T_rate = 2.883*(2018) - 5666;   
T_cost = T_hrd*T_rate;
%fprintf('\n Tooling costs = $%g ', T_cost);
%--------------------------------------------------------------------------
% Manufacturing Cost
L_hr = 7.37*(W_e^.82)*(V_max^.484)*(Q^.641); %Labor hours, same units as above
L_rate = 2.316*(2018)-4552;
L_cost = L_hr*L_rate;
%fprintf('\n Manufacturing Cost = $%g ', L_cost);
%--------------------------------------------------------------------------
% Quality Control Cost
QC_hr = 0.13*L_hr;	%Quality control hours
QC_rate = 2.6*(2018) -5112;
QC_cost = QC_hr*QC_rate;
%fprintf('\n Quality Control Cost = $%g ', QC_cost);
%--------------------------------------------------------------------------
% Material Costs
M_cost = (16.39*(W_e^0.921)*(V_max^.621)*(Q^0.799))*CPI_equiv;
%fprintf('\n Material Costs = $%g ', M_cost);
%--------------------------------------------------------------------------
% Propulsion & Avionics Cost
P_per = ne*2306*(0.043*T_sl + 243.3*M_max + 0.969*Temp_Inlet - 2228)*CPI_equiv; %TSL is max sea level thrust in lbs, Temp is in Rankine
P_cost = Qd*P_per;
%fprintf('\n Propulsion & Avionics Cost = $%g ', P_cost);
%--------------------------------------------------------------------------
RTDE_Cost = AE_cost + DSC + FTC + T_cost + L_cost + QC_cost + M_cost + P_cost;
%fprintf('\n Total RTDE Cost = $%g ', RTDE_Cost);
%--------------------------------------------------% 

%% ========================================================================
% DOC:
%--------------------------------------------------------------------------
Tto =  T_SL_max/1000; % [kN] ???
[ DOC_Cost ] = DOC_SSBJ(MTOW, W_empty, W_fuel, V_cr, Tto, ne, Range, PAX);

end
