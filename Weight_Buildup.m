function [weights, weight_fractions] = Weight_Buildup( num_pass, num_crew, V_cr, M_cr, range, SFC, LD_cr)
% Purpose: parametrically determine the weight build-up for the vehicle using weight fraction analysis
% ALL WEIGHTS IN LBF. CRUISE VELOCITY IN MPH. RANGE IN MILES.
%--------------------------------------------------------------------------
% Fixed weight calculation
W_pass = 205 * num_pass; % lbf, based on Sadraey p. 97
W_lug = 140 * num_pass;
W_pl = W_pass + W_lug;
W_crew = num_crew * 200;
W_fixed = W_pl + W_crew;
%--------------------------------------------------------------------------
% Weight fractions
WF_to = 0.98;                                           % Weight fraction for takeoff/taxi (empirical)
if M_cr > 1
    WF_climb = 0.991 - 0.007 * M_cr - 0.01*(M_cr)^2;    % Weight fraction if climbing from Mach .1 to Mach 1.4 (Nicolai)
    WF_accel = WF_climb / WF_climb;                     % Weight fraction for the accel from Mach .9 to Mach 1.4
else
    WF_climb = 1.0065 - .0325 * M_cr;                   % Weight fraction from climb to Mach 0.9
    WF_accel = 1;
end
WF_cruise = exp((-range*SFC)/(V_cr*LD_cr));             % Weight fraction for cruise segment
WF_des = 0.99;                                          % Weight fraction for descent (empirical)
WF_land = 0.997;                                        % Weight fraction for landing/taxi back (empirical)
%--------------------------------------------------------------------------
% Landing to takeoff weight ratio
WF_total = WF_to*WF_climb*WF_accel*WF_cruise*WF_des*WF_land;
%--------------------------------------------------------------------------
% Fuel weight to takeoff weight ratio from the calculated weight ratios with reserve fuel
Wf_Wto = 1.05*(1-WF_total);
%--------------------------------------------------------------------------
% Empty weight to MTOW ratio, from Plumley (range must be in between 0.3 and 0.45
We_Wto = 0.3;
if Wf_Wto + We_Wto > 1
    % If the fuel weight fraction + empty weight fraction is more than 1, the script won't work 
    err_msg = sprintf('The fuel weight fraction + empty weight fraction is more than 1.0 \n Modify the empty weight fraction!');
    error('sytnax:requirements','\n\n %s \n\n',err_msg);
end
%--------------------------------------------------------------------------
W_to = W_fixed/(1 - Wf_Wto - We_Wto);
%--------------------------------------------------------------------------
W_fuel = Wf_Wto*W_to;                               
W_land = WF_total*W_to;                             
W_empty = W_to - W_fuel - W_fixed; % Roskam Table 2.14 Vol 1
%--------------------------------------------------------------------------
W_payload = struct('Passengers', W_pass,'Luggage', W_lug, 'Crew', W_crew);
weights = struct('W_empty', W_empty, 'W_payload', W_payload, 'W_fuel', W_fuel, 'W_to', W_to, 'W_land', W_land);
weight_fractions = struct('WF_to', WF_to, 'WF_climb', WF_climb, 'WF_accel', WF_accel, 'WF_cruise', WF_cruise, 'WF_des', WF_des, 'WF_land', WF_land, 'WF_total', WF_total);

end
