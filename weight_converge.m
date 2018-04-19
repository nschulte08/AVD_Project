%{
    Weight convergence function
%--------------------------------------------------------------------------
INPUTS:
MTOW            = max take off weight (N)
alt_cr_super    = supersonic cruise altitude (m)
M_cr_super      = supersonic cruise Mach #
range_super     = supersonic cruise range requirement
TSFC_super      = supersonic cruise TSFC (1/hr)
LD_cruise       = empirical supersonic cruise lift to drag ratio
num_crew        = number of crew
num_pass        = number of passengers
ThrustLoading   = design point thrust loading
%--------------------------------------------------------------------------
OUTPUTS:
WEIGHTS = data structure with all the weights in Newtons
%--------------------------------------------------------------------------
Last modified: 04/18/2018
% =========================================================================
%}
function [WEIGHTS] = weight_converge(MTOW, alt_cr_super, M_cr_super, range_super, TSFC_super, LD_cruise, num_crew, num_pass, ThrustLoading)
%% ========================================================================
T_W = ThrustLoading;
MTOW_e = convforce(MTOW,'N','lbf'); % first guess at max weight, lbf

[~, ~, ~, ~, son_climb_super, ~, ~, ~, ~, ~] = ATMO(alt_cr_super, 'M');
V_cr_super = M_cr_super*son_climb_super; % cruise velocity, (m/s)

[weights, wt_frac] = Weight_Buildup(MTOW_e, num_pass, num_crew, convvel(V_cr_super,'m/s','mph'), M_cr_super, convlength(range_super,'m','mi'), TSFC_super, LD_cruise);
%--------------------------------------------------------------------------
% Coleman weight convergence parameters and other constant weights:
W_sys = 20000;  % systems weight, lbf (Torenbeek, p288)
W_op = 0;       % operational items weight, lbf
E_TW = 10;      % engine thrust to engine weight ratio (Coleman, p383)
mu_a = 0;       % OEW margin

W_crw = weights.W_payload.Crew;                                     % crew weight, lbf
W_pay = weights.W_payload.Passengers + weights.W_payload.Luggage;   % payload weight, lbf (pax + baggage)
%--------------------------------------------------------------------------
W_fuel = weights.W_fuel;            % first guess at fuel weight, lbf
W_str = weights.W_empty - W_sys;    % first guess at structural weight, lbf
OEW(1) = MTOW_e  - W_fuel - W_pay;  % first guess at operating empty weight, lbf
TOGW(1) = MTOW_e;                   % first guess at take off gross weight, lbf
%% ========================================================================
% weight convergence:
tol = 1;          % tolerance (10 lbf)
delta   = tol + 1; % initilize > tol
delta_2 = tol + 1; % initilize > tol
m = 1;             % TOGW index
while delta > tol    
    %% ====================================================================
    % OEW convergence:
    n = 1; % OWE index
    while delta_2 > tol
        
        OWE = OEW(n) + W_pay + W_crw;  % operating weight empty, OEW + W_pay + W_crw, lbf
        WR = TOGW(m)/OWE;              % weight ratio, TOGW/OWE
        f_sys = W_sys/OEW(n);          % ratio of systems weight to OEW
        
        OEW(n+1) = (W_str + W_sys + W_op + ((T_W*WR)/E_TW)*(W_pay + W_crw))/(1/(1 + mu_a) - f_sys - (T_W*WR)/E_TW);
        
        delta_2 = abs(OEW(n+1) - OEW(n));
        n = n+1;
        
    end
    
    OEW = OEW(n);
    %% ====================================================================
	TOGW(m+1) = (OEW + W_pay)/(1 - wt_frac.WF_total);
    W_fuel = TOGW(m+1)*wt_frac.WF_total; % fuel weight, lbf
    W_str = OEW - W_sys;  % structural weight, lbf
    
    [~, wt_frac] = Weight_Buildup(TOGW(m+1), num_pass, num_crew, convvel(V_cr_super,'m/s','mph'), M_cr_super, convlength(range_super,'m','mi'), TSFC_super, LD_cruise);
    
    delta = abs(TOGW(m+1) - TOGW(m));
    m = m+1;
    
end
MTOW = convforce(TOGW(m),'lbf','N');
W_fuel = convforce(W_fuel,'lbf','N');
% final weight fractions:
[~, wt_frac] = Weight_Buildup(TOGW, num_pass, num_crew, convvel(V_cr_super,'m/s','mph'), M_cr_super, convlength(range_super,'m','mi'), TSFC_super, LD_cruise);

%% ========================================================================
% final weights:
WEIGHTS.MTOW = MTOW;
WEIGHTS.W_fuel = W_fuel;
WEIGHTS.W_empty = convforce(OEW,'lbf','N'); % (N)
WEIGHTS.W_payload_total = convforce(W_pay,'lbf','N'); % (N)
%--------------------------------------------------------------------------
WEIGHTS.W_to_end             = MTOW*wt_frac.WF_to;                                                                      % Wt at end of TO, start of climb (N)
WEIGHTS.W_climb_end_super    = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel;                                    % Wt at end of climb, start of cruise (N)
WEIGHTS.W_climb_avg_super    = 0.5*(WEIGHTS.W_to_end + WEIGHTS.W_climb_end_super);                                      % average climb weight to supersonic altitude (N)
WEIGHTS.W_cruise_start_super = WEIGHTS.W_climb_end_super;                                                               % Wt at beginning of cruise (N)
WEIGHTS.W_cruise_end_super   = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise;                  % Wt at end of cruise (N)
WEIGHTS.W_cruise_avg_super   = 0.5*(WEIGHTS.W_climb_end_super + WEIGHTS.W_cruise_end_super);                            % Average cruise wt (N)
WEIGHTS.W_descend_end_super  = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise*wt_frac.WF_des;   % Wt at end of descent (N)
WEIGHTS.W_descend_avg_super  = 0.5*(WEIGHTS.W_cruise_end_super + WEIGHTS.W_descend_end_super);                          % Average descent wt (N)
WEIGHTS.W_climb_end_sub      = MTOW*wt_frac.WF_to*wt_frac.WF_climb;                                                     % Wt at end of climb, start of cruise (N)
WEIGHTS.W_climb_avg_sub      = 0.5*(WEIGHTS.W_to_end + WEIGHTS.W_climb_end_sub);                                        % average climb weight to subsonic altitude (N)
WEIGHTS.W_cruise_start_sub   = WEIGHTS.W_climb_end_sub;                                                                 % Wt at beginning of cruise (N)
WEIGHTS.W_cruise_end_sub     = WEIGHTS.W_climb_end_sub - W_fuel;                                                        % Wt at end of cruise (N) [approximate]
WEIGHTS.W_cruise_avg_sub     = 0.5*(WEIGHTS.W_climb_end_sub + WEIGHTS.W_cruise_end_sub);                                % Average cruise wt (N)
WEIGHTS.W_descend_end_sub    = WEIGHTS.W_cruise_end_sub*0.99;                                                           % Wt at end of descent (N)
WEIGHTS.W_descend_avg_sub    = 0.5*(WEIGHTS.W_cruise_end_sub + WEIGHTS.W_descend_end_sub);                              % Average descent wt (N)
WEIGHTS.W_land               = MTOW*(wt_frac.WF_total/wt_frac.WF_land);                                                 % landing weight (N)
WEIGHTS.W_cruise_avg_avg     = (WEIGHTS.W_cruise_avg_sub + WEIGHTS.W_cruise_avg_super)/2;                               % used for op_envelope and Ta vs Tr plot

end
