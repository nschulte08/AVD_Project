%{
function to build the operational envelope
---------------------------------------------------------------------------
INPUTS:
W_cruise_avg    = average cruise weight [N]
T_max           = max thrust
Sref            = reference area [m^2]
SM              = static margin
b_unswept       = unswept span [m]
TR              = taper ratio
CL_max          = max lift coefficient
ne              = number of engines
---------------------------------------------------------------------------
OUTPUTS:

===========================================================================
%}
function op_envelope(W_cruise_avg, Sref, SM, b_unswept, TR, CL_max, ne)
%% ========================================================================
% Operational Envelope
%--------------------------------------------------------------------------
q_max = 86184.414; % [Pa] max. dyn. pressure (p.104 Nicolai) = 1800 psf 
%--------------------------------------------------------------------------
% Altitude and Velocity range:
h = 0:5000:10e4; % [m]
%V = 0:5:800;     % [m/s]
M = [0.1:0.1:0.9,1.1:0.1:3];
%M = 1.1:0.01:3;
%--------------------------------------------------------------------------
for n = 1:length(h)

   [~, ~, ~, rho, a, ~, ~, ~, ~, ~] = ATMO(h(n), 'M');
 
for m = 1:length(M)
     V(m) = M(m)*a;
%--------------------------------------------------------------------------
% Thrust available:
[T_A_single_engine, ~] = Propulsion(M(m), h(n));
T_A = T_A_single_engine*ne; % [N] total takeoff thrust
%--------------------------------------------------------------------------
CL(n,m) = W_cruise_avg/(0.5*rho*V(m)^2*Sref); % steady level flight
%--------------------------------------------------------------------------
% Calculate sweep:
M_perp = 0.7;                        % perp Mach #, variable to iterate(?)
sweep_deg = acosd(M_perp/M(m));      % This has to be limited to 70 ish degrees!
if sweep_deg > 70                    % Limit the sweep angle
    sweep_deg = 70;
end
b_swept = b_unswept*cosd(sweep_deg); % Span at sweep angle [m]
AR = b_swept^2/Sref;           % Swept aspect ratio
%--------------------------------------------------------------------------
[CD, CD0] = aerofunk_drag_2(h(n), M(m), Sref, CL(n,m), SM, AR, TR);
D = CD*0.5*rho*V(m)^2*Sref; % [N] Drag (sluf)

Ps(n,m) = V(m)*(T_A - D)/W_cruise_avg; % Specific Excess Power
%--------------------------------------------------------------------------
end;
 V_stall_op(n) = sqrt(2*W_cruise_avg/(CL_max*rho*Sref)); % stall boundary
 M_stall_op(n) = V_stall_op(n)/a;
 V_q_lim(n) = sqrt(2*q_max/rho); % dynamic pressure limit
 M_q_lim(n) = V_q_lim(n)/a;
end;
%--------------------------------------------------------------------------
figure('Name','Operational Envelope','NumberTitle','off');
hold on;
vals = [1,1];
%vals = 0:50;
[C1, h1] = contour(M,h,Ps,vals);
set(h1, 'LineWidth', 2.5)
set(h1, 'Color', 'k')
title('Operational Envelope','FontSize',18);
xlabel('Mach Number','FontSize',12);
ylabel('Altitude (m)','FontSize',12);
% xlim([0, 2.5]);
% ylim([0, 2e4]);
plot (M_stall_op,h,':k','LineWidth',2.5);
plot (M_q_lim,h,'--k','LineWidth',2.5);
%plot (M_sub_R,h_sub_R,'k+','LineWidth',2.5,'LineStyle','None','MarkerSize',16);
%plot (M_super_R,h_super_R,'k*','LineWidth',2.5,'LineStyle','None','MarkerSize',16);
%legend('Thrust Limit','Stall Boundary', 'Dyn. Press. Limit', 'Subsonic Cruise', 'Supersonic Cruise','Location','best');
legend('Thrust Limit','Stall Boundary', 'Dyn. Press. Limit', 'Location','best');
%--------------------------------------------------------------------------
