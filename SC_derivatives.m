%{ 
===========================================================================
Stability Derivatives function

===========================================================================
INPUTS:
h = altitude [m]
M = Mach #
AR = apect ratio
Lambda = sweep angle [degrees]
S_ref = referrence area [m^2]
c_r = Root Chord [m]
b = span [m]
TR = Taper Ratio
===========================================================================
OUTPUTS:
CMa     = longitudinal stability derivative
Cl_beta = lateral stability derivative
Cn_beta = directional stability derivative
Cn_dr   = rudder control power
S_VT    = total required vertical tail area [m^2]
l_VT    = distance (in x-direction) between cg location and 1/4 chord of VT mac [m]
===========================================================================
%}
function [CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = SC_derivatives(h, M, AR, Lambda, S_ref, b, c_r, TR)
%--------------------------------------------------------------------------
% geometry:
c_bar = ((2/3)*c_r*((1 + TR + TR^2)/(1 + TR)));  % Nicolai Pg 580
%--------------------------------------------------------------------------
% atmospheric and aerodynamic properties:
[~, ~, ~, rho_atm, son_atm, ~, ~, ~, ~, ~] = ATMO(h, 'M');
%V  = M*son_atm; % [m/s] Real Velocity
%q_bar = 0.5*rho_atm*V^2; % [Pa]
%D = CD*S_ref*q_bar; % [N]

SM = 1; % static margin ---> could do a trade study with this
[~, ~,  CL, CD, CD0, CM, CLa, CMa, sweep_rad, alpha_trim] = aerofunky3(h, M, AR, TR, S_ref, SM);

%% ========================================================================
% vertical tail sizing:
C_VT = 0.09;            % Nicolai table 11.8 (buisness aircraft)
lVT_SVT = C_VT*b*S_ref; % [m^3] product of l_VT and S_VT (Nicolai p.286)
%--------------------------------------------------------------------------
% plot solutions:
l_VT_plot = 1:0.1:5;            % [m] distance (in x-direction) between cg location and 1/4 chord of VT mac
S_VT_plot = lVT_SVT./l_VT_plot; % [m^2] VT area

figure_name = 'Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location';
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
plot(l_VT_plot,S_VT_plot, 'k', 'LineWidth',3);
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
title('Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location','FontSize',18);
grid on
%--------------------------------------------------------------------------
% inital estimate:
l_VT = (3/8)*b*cosd(Lambda); % this is a function of sweep
S_VT = lVT_SVT/l_VT;
CL_a_VT = 2*pi; % [1/rad] vertical tail lift curve slope (approximation for now)

fprintf('\n\n ============================================================= \n');
fprintf('\n Vertical tail size estimate:');
fprintf('\n\n Total required vertical tail area:');
fprintf('\n S_VT = %g [m^2]', S_VT);
fprintf('\n\n Distance between cg location and 1/4 chord of VT mac:');
fprintf('\n l_VT = %g [m]', l_VT);
fprintf('\n\n ============================================================= \n');

%% ========================================================================
% Longitudinal stability: (determined from aero function)
%--------------------------------------------------------------------------
% CM_de = ; % elevator control power -------------------> need to add this
% CM_da = ; % effect of ailerons on pitching moment ----> need to add this

fprintf('\n\n ============================================================= \n');
fprintf('\n Longitudinal Stability Derivatives:');
fprintf('\n\n CM_alpha = %g [1/rad]', CMa);
fprintf('\n CM_alpha = %g [1/deg]'  , CMa*pi/180);
fprintf('\n\n CL_alpha = %g [1/rad]', CLa);
fprintf('\n CL_alpha = %g [1/deg]'  , CLa*pi/180);
fprintf('\n\n ============================================================= \n');
%% ========================================================================
% Lateral stability:
%--------------------------------------------------------------------------
% wing contribution:
Cl_b_basic = -1/AR;  % completely fictional equation ---> need to find this ---> use DATCOM handbook?? (via Nicolai's suggestion)
Cl_b_Lambda = -0.05; % due to sweep? ---> what equation is this?? fig. 21.10...?
Cl_b_Gamma = 0;     % no dihedral
Cl_b_wing = Cl_b_basic + Cl_b_Lambda + Cl_b_Gamma;

% vertical tail contribution:
z_v = 1; % [m] vertical distance between cg location and vertical tail ac (place holder)
var_2115 = 0.724 + (3.06*(S_VT/S_ref))/(1 + cosd(Lambda)) + 0.009*AR; % Nicolai eq. 21.15
Cl_b_VT = -CL_a_VT*var_2115*(S_VT/S_ref)*(z_v/b);

Cl_beta = Cl_b_wing + Cl_b_VT; % lateral stability derivative

%Cl_da = ; % aileron control power ---> use Roskam

fprintf('\n\n ============================================================= \n');
fprintf('\n Lateral Stability Derivatives:');
fprintf('\n\n Cl_beta = %g [1/rad]', Cl_beta);
fprintf('\n Cl_beta = %g [1/deg]'  , Cl_beta*pi/180);
fprintf('\n\n ============================================================= \n');

%% ========================================================================
% Directional stability:
%--------------------------------------------------------------------------
% wing contribution:
x = 1; % [m] distance between cg location and wing ac (place holder)
Cn_b_wing = CL^2*(1/(4*pi*AR) - tand(Lambda)/(pi*AR*(AR + 4*cosd(Lambda)))*(cosd(Lambda) - AR/2 - AR^2/(8*cosd(Lambda)) + 6*x*sind(Lambda)/(c_bar*AR))); % Nicolai eq. 21.22

% vertical tail contribution:
Cn_b_VT = C_VT*CL_a_VT*var_2115; % Nicolai eq. 21.21

% directional stability derivative:
Cn_beta = Cn_b_wing + Cn_b_VT; % Nicolai eq. 21.20

tau = 0.6; % placeholder (Nicolai Figure 21.14) ---> this fig could be digitized and the data interpolated
Cn_dr = 0.9*CL_a_VT*C_VT*tau; % rudder control power (Nicolai eq. 21.26)

fprintf('\n\n ============================================================= \n');
fprintf('\n Directional Stability Derivatives:');
fprintf('\n\n total:');
fprintf('\n Cn_beta = %g [1/rad]', Cn_beta);
fprintf('\n Cn_beta = %g [1/deg]', Cn_beta*pi/180);
fprintf('\n\n vertical tail contribution:');
fprintf('\n Cn_b_VT = %g [1/rad]', Cn_b_VT);
fprintf('\n Cn_b_VT = %g [1/deg]', Cn_b_VT*pi/180);
fprintf('\n\n rudder control power:');
fprintf('\n Cn_dr = %g [1/rad]', Cn_dr);
fprintf('\n Cn_dr = %g [1/deg]', Cn_dr*pi/180);
fprintf('\n\n ============================================================= \n');

%--------------------------------------------------------------------------
dr_plot = 0:2:10;     % [deg] rudder deflection
beta_plot = -10:10; % [deg] side slip

    
figure_name = 'Directional Stability';
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
for m = 1:length(dr_plot)
    Cn_plot(:,m) = Cn_beta.*beta_plot.*pi./180 + Cn_dr*dr_plot(m)*pi/180;
end
plot(beta_plot, Cn_plot(:,1),  'k', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,2),  'b', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,3),  'c', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,4),  'r', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,5),  'm', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,6),  'g', 'LineWidth',3);

hold off
xlabel('\beta (deg)'  ,'FontSize',18);
ylabel('C_n','FontSize',18);
title('Yawing Moment vs Sideslip','FontSize',18);
legend('\delta_r = 0 [deg]','\delta_r = 2 [deg]','\delta_r = 4 [deg]','\delta_r = 6 [deg]','\delta_r = 8 [deg]','\delta_r = 10 [deg]','Location', 'best');
grid on

%% ========================================================================

end

