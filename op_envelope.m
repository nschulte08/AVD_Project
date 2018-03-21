%{
function to build the operational envelope
---------------------------------------------------------------------------
INPUTS:
W_cruise_avg    = average cruise weight [N]
T_max           = max thrust
Sref            = reference area [m^2]
SM              = static margin
AR              = aspec ratio
TR              = taper ratio
CL_max          = max lift coefficient
---------------------------------------------------------------------------
OUTPUTS:

===========================================================================
%}
function op_envelope(W_cruise_avg, T_max, Sref, SM, AR, TR, CL_max)
%% ========================================================================
% Operational Envelope
%--------------------------------------------------------------------------
q_max = 86184.414; % [Pa] max. dyn. pressure (p.104 Nicolai) = 1800 psf 
%--------------------------------------------------------------------------
% Altitude and Velocity range:
h = 0:500:25000; % [m]
%V = 0:5:800;     % [m/s]
M = [0:0.01:0.9,1.1:0.01:3];
%--------------------------------------------------------------------------
for n = 1:length(h)

   [~, ~, ~, rho, a, ~, ~, ~, ~, ~] = ATMO(h(n), 'M');
 
for m = 1:length(M)
     V(m) = M(m)*a;
%--------------------------------------------------------------------------
% Thrust available:
%{
T_A_data = csvread('Thrust_NEW.csv'); % [lbf] per engine
T_A_data = T_A_data*4.44822*3;        % [N] all three engines
MM = 0:0.1:2;
hh = 0:1000:100000; % [ft]
hh = hh*0.3048;     % [m]

Thrust_A = interp2(MM,hh,T_A_data, M(m), h(n)); % [N] thrust available
%}
T_A = Thrust(h(n), M(m), T_max);
%--------------------------------------------------------------------------
CL(n,m) = W_cruise_avg/(0.5*rho*V(m)^2*Sref); % steady level flight
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
legend('Thrust Limit','Stall Boundary', 'Dyn. Press. Limit', 'Subsonic Cruise', 'Supersonic Cruise','Location','best');
%--------------------------------------------------------------------------
