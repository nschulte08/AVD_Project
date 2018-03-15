function [S_descend, dt_descend] = perf_descent(alt_sub_cr, M_cr)
%This function returns the performance parameters for the OFW during
%descent

%Inputs:
%alt_sub_cr:    Altitude of subsonic cruise
%M_cr:          Cruise Mach number

%Outputs:
%S_descend:     Descent Distance Covered
%dt_descend:    Time of descent phase

%Local Inputs:
M_sub = 0.95;
[~,a,~,~] = atmosisa(alt_sub_cr);
V_sub = a * M_sub;


%----------------------------  Descent  ---------------------------------
V_descend = V_sub;                                  % [m/s] descent velocity
gamma_descend = -atand(1000/(3*5280));                % 3:1 rule
ROD = V_descend*sind(gamma_descend);                  % [m/s] rate of descent
%--------------------------------------------------------------------------
if M_cr < 1
    S_descend = h_sub_R/abs(tand(gamma_descend));     % [m] range covered during descent
    dt_descend = S_descend/V_descend;             % [s] time to descend
    
    fprintf('\n\n ============================== Descent Results  ============================== \n');
    fprintf('\n Descent from Subsonic Cruising Altitude: (h = %g [m])', h_sub_R);
    fprintf('\n Descent Angle:         gamma = %g [deg] ', gamma_descend);
    fprintf('\n Rate of Descent:         ROD = %g [m/s] = %g [ft/s]', ROD, ROD*3.2808399);
    fprintf('\n Descent Velocity:          V = %g [m/s] = %g [ft/s]', V_descend, V_descend*3.2808399);
    fprintf('\n\n -------------------------------------------------------------------- ');
    fprintf('\n Range Covered During Descent:   R = %g [km] = %g [miles]', S_descend/1000, S_descend*0.000621371);
    fprintf('\n Time to Descend:               dt = %g [min]', dt_descend/60);
    fprintf('\n\n -------------------------------------------------------------------- ');

else
    S_descend = h_super_R/abs(tand(gamma_descend)); % [m] range covered during descent
    dt_descend = S_descend/V_descend;         % [s] time to descend

    fprintf('\n --------------------------------------------------------------------\n ');
    fprintf('\n Descent from Supersonic Cruising Altitude: (h = %g [m])', h_super_R);
    fprintf('\n Descent Angle:         gamma = %g [deg] ', gamma_descend);
    fprintf('\n Rate of Descent:         ROD = %g [m/s] = %g [ft/s]', ROD, ROD*3.2808399);
    fprintf('\n Descent Velocity:          V = %g [m/s] = %g [ft/s]', V_descend, V_descend*3.2808399);
    fprintf('\n\n -------------------------------------------------------------------- ');
    fprintf('\n Range Covered During Descent:   R = %g [km] = %g [miles]', S_descend/1000, S_descend*0.000621371);
    fprintf('\n Time to Descend:               dt = %g [min]', dt_descend/60);
    fprintf('\n\n ============================================================================== \n');
end

end

