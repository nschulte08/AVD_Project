%{
    Calculate the RDTE and DOC costs
---------------------------------------------------------------------------
INPUTS:
W_e     =
V_max   =
T_sl    =
M_max   =
W_to    =
W_f     =
Range   = total range [miles???????????]
t_cl    = Climb time [hr]
t_cr    = Cruise time [hr]
t_de    = Descent time [hr]
W_A     =
N_e     = number of engines???????????
---------------------------------------------------------------------------
OUTPUTS:
RTDE_Cost = research developement and E cost [$]
---------------------------------------------------------------------------
%}
function [ RTDE_Cost ] = costfunky( W_e, V_max, T_sl, M_max, W_TO, W_f, Range, t_cl, t_de, t_cr, W_A, N_e)

fprintf('\n\n ================================= Cost Results ================================= \n');

%% RDTE
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
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cost of airframe engineering = $%g ', AE_cost);
%--------------------------------------------------------------------------
% Development Support Cost (D)
DSC = (66*(W_e^0.63)*(V_max^1.3))*CPI_equiv;        %Cost for development support, We in lbs, Vmax in knots
fprintf('\n Cost for development support = $%g ', DSC);
%--------------------------------------------------------------------------
% Flight Test Op Cost (F)
FTC = (1852*(W_e^0.325)*(V_max^0.822)*(Qd^1.21))*CPI_equiv;	%Cost for flight test, same units as above
fprintf('\n Cost for flight test = $%g ', FTC);
%--------------------------------------------------------------------------
%Tooling costs
T_hrd = 5.99*(W_e^.777)*(V_max^.696)*(Q^.263);	%Tooling hours, same units as above
T_rate = 2.883*(2018) - 5666;   
T_cost = T_hrd*T_rate;
fprintf('\n Tooling costs = $%g ', T_cost);
%--------------------------------------------------------------------------
% Manufacturing Cost
L_hr = 7.37*(W_e^.82)*(V_max^.484)*(Q^.641); %Labor hours, same units as above
L_rate = 2.316*(2018)-4552;
L_cost = L_hr*L_rate;
fprintf('\n Manufacturing Cost = $%g ', L_cost);
%--------------------------------------------------------------------------
% Quality Control Cost
QC_hr = 0.13*L_hr;	%Quality control hours
QC_rate = 2.6*(2018) -5112;
QC_cost = QC_hr*QC_rate;
fprintf('\n Quality Control Cost = $%g ', QC_cost);
%--------------------------------------------------------------------------
% Material Costs
M_cost = (16.39*(W_e^.921)*(V_max^.621)*(Q^.799))*CPI_equiv;
fprintf('\n Material Costs = $%g ', M_cost);
%--------------------------------------------------------------------------
% Propulsion & Avionics Cost
P_per = N_e*2306*(0.043*T_sl + 243.3*M_max +.969*Temp_Inlet - 2228)*CPI_equiv; %TSL is max thrust in lbs, Temp is in Rankine
P_cost = Qd*P_per;
fprintf('\n Propulsion & Avionics Cost = $%g ', P_cost);
%--------------------------------------------------------------------------
RTDE_Cost = AE_cost + DSC + FTC + T_cost + L_cost + QC_cost + M_cost + P_cost;
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Total RTDE Cost = $%g ', RTDE_Cost);
fprintf('\n -------------------------------------------------------------------------------- ');
%--------------------------------------------------% 
%% Direct Operating Cost
% Inputs
U_ann = [500 1500];                 % Annual utilization (long haul) [hr] {1 - p77}
t_gm = 0.51e-6*W_TO + 0.125;        % Time spent during ground maneuvers [hr] {1 - p71}
t_bl = t_gm + t_cl + t_cr + t_de;   % Total block flight time [hr] {1 - p71}
V_bl = Range/t_bl;                  % Total block flight velocity [nmi/h] {1 - 71}
%-----------------------------------------------------------------------------------------------%
% Crew Cost
% Inputs
K1 = 0.26;              % Captain's benefits/insurance factor {1}
K2 = 0.26;              % Co-pilot's benefits/insurance factor {1}
SAL1 = [30000 72000];   % PLACEHOLDER Captain's annual salary [USD]
SAL2 = [22000 52000];   % PLACEHOLDER Co-pilot annual salary [USD]
AH1 = 750;              % PLACEHOLDER Number of captain's flight hours per year {1}
AH2 = 750;              % PLACEHOLDER Number of co-pilot's flight hours per year {1}
TEF1 = 11;              % PLACEHOLDER Captain travel expense factor [USD/hr] {1}
TEF2 = 11;              % PLACEHOLDER Co-pilot travel expense factor [USD/hr] {1}

% Crew cost (Assuming 2 crew) [USD/nmi] {1}:
C_crew = ((1 + K1)/V_bl)*(SAL1/AH1) + TEF1/V_bl + ((1 + K2)/V_bl)*(SAL2/AH2) + TEF2/V_bl;
%-----------------------------------------------------------------------------------------------%
% Petroleum, Oil, and Lubricants Cost
% Inputs
FP = 5.20;      % PLACEHOLDER Fuel price (Jet A) [USD/ga]
FD = 6.74;      % PLACEHOLDER Fuel density (Jet A) [lbm/ga] {1}
% Fuel and oil cost, assuming oil is 5% of overall POL costs [USD/nmi] {1}:
C_pol = 1.05*(W_f/Range)*(FP/FD);
%-----------------------------------------------------------------------------------------------%
% Airframe Insurance Cost
% Inputs
f_ins_hull = [0.005 0.03];  % PLACEHOLDER Annual hull insurance rate [USD/aircraft price/aircraft/year] {1}
AMP = 120e6;                % PLACEHOLDER Aircraft market price [USD]
% Analysis
C_ins = f_ins_hull*AMP/(U_ann*V_bl); % Airframe insurance cost [USD/nmi] {1}
%-----------------------------------------------------------------------------------------------%
% Direct Operating Cost of Flying
DOC_fly_min = min(C_crew) + min(C_pol) + min(C_ins);    % Minimum direct operating cost of flying [USD/nmi]
DOC_fly_max = max(C_crew) + max(C_pol) + max(C_ins);    % Maximum direct operating cost of flying [USD/nmi]
DOC_fly = [DOC_fly_min DOC_fly_max];
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n % Minimum direct operating cost of flying  = %g [USD/nmi]', DOC_fly_min);
fprintf('\n % Maximum direct operating cost of flying  = %g [USD/nmi]', DOC_fly_max);
%-----------------------------------------------------------------------------------------------%
% Direct Operating Cost of Maintenance:

% Airframe and Non-Engine-Systems Maintenance
 % Inputs
 MHR_flt_lab = 3.0;         % PLACEHOLDER Airframe and systems mainteance man hours per flight hour {1}
 R_l_ap = [5.00 18.40];     % PLACEHOLDER Aircraft maintenance labor rate [USD/hr] {1}
 % Analysis
 t_flt = t_cl + t_cr + t_de;                % Average flight time [hr]
 MHR_bl_lab = MHR_flt_lab*(t_flt/t_bl);     % Airframe and systems mainteance man hours per block hour {1}
 C_lab_ap = 1.03*MHR_bl_lab*R_l_ap/V_bl;    % Airframe and non-engine-systems maintenance labor costs [USD/nmi] {1}
%--------------------------------------------------------------------------
% Engine maintenance
 % Analysis
 MHR_bl_eng = 3.0*0.067*W_A/1000;                   % Engine maintenance hours per block hour per engine {1}
 R_l_eng = [5.00, 18.40];                           % Engine maintenance labor rate [USD/hr] {1}
 C_lab_eng = 1.03*1.3*N_e*MHR_bl_eng*R_l_eng/V_bl;  % Engine maintenance labor costs [USD/nmi] {1}
%--------------------------------------------------------------------------
% Airframe and Non-Engine-Systems Materials Maintenance
 % Inputs
 C_mat_ap_bl = 300;     % PLACEHOLDER Airframe and non-engine-systems maintenance material cost per block hour [USD/hr] {1}
 % Analysis
 C_mat_ap = 1.03*C_mat_ap_bl/V_bl;  % Airframe and non-engine-systems materials maintenance costs [USD/nmi] {1}
 %--------------------------------------------------------------------------
% Engine Materials Maintenance
 % Inputs
 EP = 1.5e6;            % PLACEHOLDER Engine price [USD]
 ESPPF = 1.5;           % PLACEHOLDER Engine spare parts price factor {1}
 H_em = [3000 5000];    % Period between engine overhauls [hr] {1}
 % Analysis
 K_H_em = 0.021*(H_em/100) + 0.769;                 % Period between engine overhaul factor {1}
 C_mat_eng_bl = (5.43e-5*EP*ESPPF - 0.47)./K_H_em;  % Engine maintenance material cost per block hour per engine [USD/hr] {1}
 C_mat_eng = 1.03*1.3*N_e*C_mat_eng_bl/V_bl;        % Engine materials maintenance costs [USD/nmi] {1}
%-------------------------------------------------------------------------- 
% Applied Maintenance Burden
 % Inputs
 f_amb_lab = [0.9, 1];     % Overhead distribution factor for labor cost {1}
 f_amb_mat = [0.3, 0.4];   % Overhead distribution factor for material cost {1}
 % Applied maintenance burden [USD/nmi] {1}:
 C_amb = 1.03*(f_amb_lab.*(MHR_bl_lab*R_l_ap + N_e*MHR_bl_eng*R_l_eng) + f_amb_mat.*(C_mat_ap_bl + N_e*C_mat_eng_bl))/V_bl;
%--------------------------------------------------------------------------
% Final Maintenance DOC Analysis
% Minimum direct operating cost of maintenance [USD/nmi] {1}:
DOC_maint_min = min(C_lab_ap) + min(C_lab_eng) + min(C_mat_ap) + min(C_mat_eng) + min(C_amb);
% Maximum direct operating cost of maintenance [USD/nmi] {1}:
DOC_maint_max = max(C_lab_ap) + max(C_lab_eng) + max(C_mat_ap) + max(C_mat_eng) + max(C_amb);
% Direct Operating Cost of Depreciation:
DOC_maint = [DOC_maint_min, DOC_maint_max];
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Minimum direct operating cost of maintenance  = %g [USD/nmi]', DOC_maint_min);
fprintf('\n Maximum direct operating cost of maintenance  = %g [USD/nmi]', DOC_maint_max);
%-----------------------------------------------------------------------------------------------%
% Airframe Depreciation Cost
 % Inputs
 F_dap = 0.85;              % PLACEHOLDER Airframe depreciation factor {1}
 ASP = [.05*AMP .4*AMP];    % PLACEHOLDER Avionics systems price per aircraft [USD] {1}
 DP_ap = 20;                % PLACEHOLDER Aircraft depreciation period [yr] {1 - p107}
 % Analysis
 C_dap = F_dap*(AMP - N_e*EP - ASP)./(DP_ap*U_ann*V_bl);    % Cost of airframe depreciation [USD/nmi] {1}
%--------------------------------------------------------------------------
% Engine Depreciation Cost
 % Inputs
 F_deng = 0.85;     % PLACEHOLDER Engine depreciation factor {1 - p107}
 DP_eng = 7;        % PLACEHOLDER Engine depreciation period {1 - p107}
 % Analysis
 C_deng = (F_deng*N_e*EP)./(DP_eng*U_ann*V_bl);     % Cost of engine depreciation [USD/nmi] {1}
%--------------------------------------------------------------------------
% Avionics Depreciation Cost
 % Inputs
 F_dav = 1;     % PLACEHOLDER Avionics system depreciation factor {1 - p107}
 DP_av = 5;     % PLACEHOLDER Avionics system depreciation period {1 - p107}
 % Analysis
C_dav = (F_dav*ASP)./(DP_av*U_ann*V_bl);   % Cost of avionics depreciation [USD/nmi] {1}
%--------------------------------------------------------------------------
% Space Parts Depreciation
 % Inputs
 F_dapsp = 0.85;    % PLACEHOLDER Aircraft spare parts depreciation factor {1 - p107}
 F_apsp = .1;       % PLACEHOLDER Aircraft spare parts factor {1 - p105}
 DP_apsp = 10;      % PLACEHOLDER Aircraft spare parts depreciation period {1 - p107}
 % Cost of spare parts depreciation [USD/nmi] {1 - p105}:
 C_dapsp = (F_dapsp*F_apsp*(AMP - N_e*EP))./(DP_apsp*U_ann*V_bl);   
%--------------------------------------------------------------------------
% Engine Spare Parts Depreciation
 % Inputs
 F_dengsp = 0.85;   % Engine spare parts depreciation factor {1 - p107}
 F_engsp = 0.5;     % Engine spare parts factor {1 - p106}
 DP_engsp = 7;      % Engine spare parts depreciation period {1 - p107}
 % Cost of engine spare parts depreciation [USD/nmi] {1 - p106}:
 C_dengsp = (F_dengsp*F_engsp*N_e*EP*ESPPF)./(DP_engsp*U_ann*V_bl);
%--------------------------------------------------------------------------
% Final Depreciation DOC Analysis 
 % Minimum direct operating cost of depreciation [USD/nmi]:
DOC_depr_min = min(C_dap) + min(C_deng) + min(C_dav) + min(C_dapsp) + min(C_dengsp);  
 % Maximum direct operating cost of depreciation [USD/nmi]:
DOC_depr_max = max(C_dap) + max(C_deng) + max(C_dav) + max(C_dapsp) + max(C_dengsp);
% Preliminary DOC estimate for DOC_lnr and DOC_fin calculation:
DOC_depr = [DOC_depr_min, DOC_depr_max];
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Minimum direct operating cost of depreciation  = %g [USD/nmi]', DOC_depr_min);
fprintf('\n Maximum direct operating cost of depreciation  = %g [USD/nmi]', DOC_depr_max);
%-----------------------------------------------------------------------------------------------%
DOC_dom0 = DOC_fly + DOC_maint + DOC_depr;  % Domestic Direct operating cost [USD/nmi] {1 - p80}
DOC_int0 = DOC_fly + DOC_maint + DOC_depr;  % International Direct operating cost [USD/nmi] {1 - p80}
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Minimum Domestic Direct operating cost  = %g [USD/nmi]', min(DOC_dom0));
fprintf('\n Maximum Domestic Direct operating cost  = %g [USD/nmi]', max(DOC_dom0));
fprintf('\n Minimum International Direct operating cost [USD/nmi]  = %g [USD/nmi]', min(DOC_int0));
fprintf('\n Maximum International Direct operating cost [USD/nmi]  = %g [USD/nmi]', max(DOC_int0));
%-----------------------------------------------------------------------------------------------%
% Direct Operating Cost of Landing, Navigation, and Registry Fees:

% Landing Fees
C_aplf = 0.002*W_TO;        % Aircraft landing fee per landing [USD/lbf] {1 - p108}
C_lf = C_aplf/(V_bl*t_bl);  % DOC due to landing fees [USD/nmi] {1 - p107}
%--------------------------------------------------------------------------
% Navigation Fees
 % Inputs
 C_apnf_dom = 0;    % Domestic (USA) flight navigation fee per aircraft per flight [USD] {1 - p109}
 C_apnf_int = 10;   % International flight navigation fee per aircraft per flight [USD] {1 - p109}
 % Analysis
 C_nf_dom = C_apnf_dom/(V_bl*t_bl);     % Cost of domestic navigation fee [USD/nmi] {1 - p108}
 C_nf_int = C_apnf_int/(V_bl*t_bl);     % Cost of international navigation fee [USD/nmi] {1 - p108}
%--------------------------------------------------------------------------
 % Registry Taxes
f_rt = 0.001 + 10e-8*W_TO;  % Aircraft size factor {1 - p109}
C_rt_dom = f_rt*DOC_dom0;   % DOC of registry taxes [USD/nmi] {1 - p109} 
C_rt_int = f_rt*DOC_int0;   % DOC of registry taxes [USD/nmi] {1 - p109}
%--------------------------------------------------------------------------
% Final Landing, Navigation, and Registry Fees DOC Analysis
% Minimum domestic direct operating cost of landing fees [USD/nmi] {1 - p107}:
DOC_lnr_dom_min = min(C_lf) + min(C_nf_dom) + min(C_rt_dom);
% Maximum domestic direct operating cost of landing fees [USD/nmi] {1 - p107}:
DOC_lnr_dom_max = max(C_lf) + max(C_nf_dom) + max(C_rt_dom);
DOC_lnr_dom = [DOC_lnr_dom_min, DOC_lnr_dom_max];
% Minimum international direct operating cost of landing fees [USD/nmi] {1 - p107}:
DOC_lnr_int_min = min(C_lf) + min(C_nf_int) + min(C_rt_int);
% Maximum international direct operating cost of landing fees [USD/nmi] {1 - p107}:
DOC_lnr_int_max = max(C_lf) + max(C_nf_int) + max(C_rt_int); 
DOC_lnr_int = [DOC_lnr_int_min, DOC_lnr_int_max];
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Minimum domestic direct operating cost of landing fees  = %g [USD/nmi]', DOC_lnr_dom_min);
fprintf('\n Maximum domestic direct operating cost of landing fees  = %g [USD/nmi]', DOC_lnr_dom_max);
fprintf('\n Minimum international direct operating cost of landing fees  = %g [USD/nmi]', DOC_lnr_int_min);
fprintf('\n Maximum international direct operating cost of landing fees  = %g [USD/nmi]', DOC_lnr_int_max);
%-----------------------------------------------------------------------------------------------%
% Direct Operating Cost of Financing:

% Final Financing DOC Analysis
DOC_fin_dom = 0.07*DOC_dom0;    % Direct operating cost of financing [USD/nmi] {1 - p109}
DOC_fin_int = 0.07*DOC_int0;    % Direct operating cost of financing [USD/nmi] {1 - p109}

% Domestic Direct operating cost [USD/nmi] {1 - p80}:
DOC_dom = DOC_fly + DOC_maint + DOC_depr + DOC_lnr_dom + DOC_fin_dom; 
% International Direct operating cost [USD/nmi] {1 - p80}:
DOC_int = DOC_fly + DOC_maint + DOC_depr + DOC_lnr_int + DOC_fin_int;

fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Minimum Domestic Direct operating cost of Financing:  = %g [USD/nmi]', min(DOC_dom));
fprintf('\n Maximum Domestic Direct operating cost of Financing:  = %g [USD/nmi]', max(DOC_dom));
fprintf('\n Minimum International Direct operating cost of Financing:  = %g [USD/nmi]', min(DOC_int));
fprintf('\n Maximum International Direct operating cost of Financing:  = %g [USD/nmi]', max(DOC_int));

%% Program Operating Cost

% Inputs
N_yr = 20;  % PLACEHOLDER Number of aircraft operating years {1}
R_ann = V_bl*U_ann;     % Total annual miles flown [nmi] {1}

% Analysis
C_OPS_dom = DOC_dom.*R_ann*N_yr;    % Domestic Direct Program Operating Costs [USD] {1}
C_OPS_int = DOC_int.*R_ann*N_yr;    % International Direct Program Operating Costs [USD] {1}
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Minimum Domestic Direct Program Operating Costs  = %g [USD/nmi]', min(C_OPS_dom));
fprintf('\n Maximum Domestic Direct Program Operating Costs  = %g [USD/nmi]', max(C_OPS_dom));
fprintf('\n Minimum International Direct Program Operating Costs  = %g [USD/nmi]', min(C_OPS_int));
fprintf('\n Maximum International Direct Program Operating Costs  = %g [USD/nmi]', max(C_OPS_int));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n\n ================================================================================ \n');

end
