%{
thrust
---------------------------------------------------------------------------
INPUTS:
cruise          = cruise points = [M_sub,h_sub, M_super,h_super];
W_cruise_avg    = average cruise weight [N]
T_max           = max thrust
Sref            = reference area [m^2]
SM              = static margin
b_unswept       = unswept span [m]
TR              = taper ratio
CL_max          = max lift coefficient
ne              = number of engines
M_perp          = perpendicular wing Mach number
---------------------------------------------------------------------------
===========================================================================
%}
function Thrust_required_and_available(cruise, W_cruise_avg, Sref, SM, b_unswept, TR, ne, M_perp)
%% ========================================================================
M_sub   = cruise(1); 
h_sub   = cruise(2); % [m]
M_super = cruise(3); 
h_super = cruise(4); % [m]

%% ========================================================================
% Operational Envelope
%--------------------------------------------------------------------------
% Altitude and Velocity range:
h = [h_sub, h_super]; % [m]
%V = 0:5:800;     % [m/s]
M = [0.1:0.1:0.9,1.1:0.1:3];
%--------------------------------------------------------------------------
for n = 1:length(h)

   [~, ~, ~, rho, a, ~, ~, ~, ~, ~] = ATMO(h(n), 'M');
 
for m = 1:length(M)
     V(m) = M(m)*a;
%--------------------------------------------------------------------------
% Thrust available:
[T_A_single_engine, ~] = Propulsion(M(m), h(n));
T_A(n,m) = T_A_single_engine*ne; % [N] total takeoff thrust
%--------------------------------------------------------------------------
CL(n,m) = W_cruise_avg/(0.5*rho*V(m)^2*Sref); % steady level flight
%--------------------------------------------------------------------------
% Calculate sweep:
sweep_deg = acosd(M_perp/M(m));      % This has to be limited to 70 ish degrees!
if sweep_deg > 70                    % Limit the sweep angle
    sweep_deg = 70;
end
b_swept = b_unswept*cosd(sweep_deg); % Span at sweep angle [m]
AR = b_swept^2/Sref;           % Swept aspect ratio
%--------------------------------------------------------------------------
[CD, CD0] = aerodynamic_drag(h(n), M(m), Sref, CL(n,m), SM, AR, TR, sweep_deg);
D = CD*0.5*rho*V(m)^2*Sref; % [N] Drag (sluf)

T_R(n,m) = D; % Thrust required
%--------------------------------------------------------------------------
end;
end;
%--------------------------------------------------------------------------
figure('Name','Thrust required and available','NumberTitle','off');
hold on;
title('Thrust required and available','FontSize',18);
xlabel('Mach Number','FontSize',12);
ylabel('Thrust (kN)','FontSize',12);
% xlim([0, 2.5]);
% ylim([0, 2e4]);
plot (M,T_R(1,:)/1000,'k','LineWidth',2.5);
plot (M,T_A(1,:)/1000,'--k','LineWidth',2.5);
plot (M,T_R(2,:)/1000,'b','LineWidth',2.5);
plot (M,T_A(2,:)/1000,'--b','LineWidth',2.5);
first = sprintf('Thrust required at %g m', h_sub);
second = sprintf('Thrust available at %g m', h_sub);
third = sprintf('Thrust required at %g m', h_super);
fourth = sprintf('Thrust available at %g m', h_super);
legend(first,second,third,fourth,'Location','best');
%--------------------------------------------------------------------------
