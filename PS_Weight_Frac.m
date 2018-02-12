clear all
%% This script is to parametrically determine the best weight & altitude 
%  for initial cruise condition

%Fuel Weight:
%inputs needed for this section: SFC, Vcruise, L/D cruise, Range

%Create altitude vector for iteration
alt = 18000;                            %Cruise altitude in meters
%Create Mach vector for iteration
Mcr = 1.2:.1:2.0;                       %Mach for cruise
%This needs to be a M,L/D vector, from Figure 4 in report    
LD_cr = [15.946, 15.1, 14.3, 13.5, 12.544, 11.75, 10.9, 10.1, 9.1];
%Create a tsfc vector based on Mach number
SFC = [0.77185, 0.8019, 0.82995, 0.8569, 0.88445, 0.9128, 0.9422, 0.9597, 0.9597];

for ii = 1:length(Mcr)
        [Tcr, acr, Pcr, rho_cr] = atmosisa(alt);        %Properties at cruise altitude
        Vcr_si = acr*Mcr(ii);                               %Cruise Velocity in m/s
        Vcr_mph = Vcr_si*2.23694;                           %Velocity in mph for formula
        
        if Mcr(ii) < 1.0
            R = 6099;                                               %Subsonic range in miles
            WF_climb = 1.0065 - .0325*(Mcr(ii));                    %Weight Fraction to Climb to Mach 0.9
            WF_accel=1.0;
        else
            R = 5466.2;                                             %Supersonic range in miles
            WF_arb = 1.0065 - .0325*.9;                             %Weight Fraction to Climb to Mach 0.9
            WF_climb = 0.991 - 0.007*Mcr(ii) - 0.01*(Mcr(ii))^2;    %Weight fraction if climbing from Mach .1 to Mach 1.4
            WF_accel = WF_climb/WF_arb;                             %This gives the wt fraction for the accel from Mach .9 to Mach 1.4   
        end
    
        WF_to = 0.98;                                       %Weight Fraction for Takeoff/Taxi (empirical)
        WF_cruise = exp((-R*SFC(ii))/(Vcr_mph*LD_cr(ii)));	%Weight fraction for cruise segment
        WF_des = 0.99;                                  	%Wt Fraction for descent (empirical)
        WF_land = 0.997;                                    %Wt Fraction for landing/taxi back (empirical)

        %Now calculate the landing to takeoff weight ratio
        WF_total = WF_to*WF_climb*WF_accel*WF_cruise*WF_des*WF_land;
 
        Wf_Wto(ii) = 1.05*(1-WF_total);        %This gives the fuelwt/takeoffwt ratio from the calculated weight ratios
end

% 
% %Find minimum of WF/WTO matrix
% Wf_min = min(Wf_Wto); 
% 
% [row,col] = find(Wf_Wto==Wf_min);
% 
% h_yay = alt(row);
% M_yay = Mcr(col);
% %%======================================================================================%%
