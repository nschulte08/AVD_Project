%{
===========================================================================
Stability Derivatives function
===========================================================================
INPUTS:
h       = altitude [m]
M       = Mach number
b       = effective span [m]
S_ref   = referrence area [m^2]
AR      = effective aspect ratio
Lambda  = sweep angle [deg]
TR      = taper ratio
CL      = lift coefficient
SM      = static margin
flight_phase = string for plot naming
===========================================================================
OUTPUTS:
CMa     = longitudinal stability derivative [1/rad]
Cl_beta = lateral stability derivative [1/rad]
Cn_beta = directional stability derivative [1/rad]
CM_de   = elevator control power [1/rad]
Cl_da   = aileron control power [1/rad]
Cn_dr   = rudder control power [1/rad]
S_VT    = total required vertical tail area [m^2]
l_VT    = distance (in x-direction) between cg location and 1/4 chord of VT mac [m]
===========================================================================
%}
function [CMa, Cl_beta, Cn_beta, CM_de, Cl_da, Cn_dr, S_VT, l_VT, vert_tail_plot] = stability(M, AR, Lambda, S_ref, b, TR, CL, SM, flight_phase)
%--------------------------------------------------------------------------
% geometry:
tmax = 2.3;         % Maximum Thickness
tcmax = 0.16;       % T/c max; variable to iterate (this is internal to aero team... ow one else uses it)
c_r = tmax/tcmax;   % [m] Root chord = max thickness / tcmax ratio
c_bar = ((2/3)*c_r*((1+TR+TR^2)/(1+TR)));  % Nicolai Pg 580

%% ========================================================================
% vertical tail sizing:
C_VT = 0.09;            % Nicolai table 11.8 (buisness aircraft)
lVT_SVT = C_VT*b*S_ref; % [m^3] product of l_VT and S_VT (Nicolai p.286)
%--------------------------------------------------------------------------
% plot solutions:
l_VT_plot = 1:0.1:b;           % [m] distance (in x-direction) between cg location and 1/4 chord of VT mac
S_VT_plot = lVT_SVT./l_VT_plot; % [m^2] VT area

vert_tail_plot = [l_VT_plot; S_VT_plot];
%{
figure_name = sprintf('Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location, for %s', flight_phase);
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
plot(l_VT_plot,S_VT_plot, 'k', 'LineWidth',3);
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
xlim([1,b]);
title_string = sprintf('%s Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location',flight_phase);
title(title_string,'FontSize',18);
grid on
fig_save('Figures', figure_name)
%}
%--------------------------------------------------------------------------
% inital estimate:
l_VT = (3/8)*b*cosd(Lambda); % this is a function of sweep
S_VT = lVT_SVT/l_VT;
CL_a_VT = 2*pi; % [1/rad] vertical tail lift curve slope (approximation for now)
%{
fprintf('\n\n ============================================================= \n');
fprintf('\n %s Vertical tail size estimate:', flight_phase);
fprintf('\n\n Total required vertical tail area:');
fprintf('\n S_VT = %g [m^2]', S_VT);
fprintf('\n\n Distance between cg location and 1/4 chord of VT mac:');
fprintf('\n l_VT = %g [m]', l_VT);
fprintf('\n\n ============================================================= \n');
%}
%% ========================================================================
% Longitudinal stability:
%--------------------------------------------------------------------------
%CMa = -1.2;
tcmax = 0.16;
Cla = 1.8*pi*(1 + 0.8*tcmax); %Airfoil Cl_alpha, Saedray pg 179  
Beta_sub = sqrt(abs(1-M^2));
nu = Cla/(2*pi/Beta_sub);
M_perp = 0.7; % perpendicular Mach #
sweep_deg = acosd(M_perp/M);
sweep_rad = sweep_deg*pi/180; % Convert to radians for formulas
CLa_sub = 2*pi*AR/(2+sqrt(4+((AR^2*Beta_sub^2)/nu^2)*(1+tan(sweep_rad)^2/Beta_sub^2))); 
CLa_super = 4/(sqrt(M^2-1))*(1-1/(2*AR*sqrt(M^2-1))); %Straight,tapered wings in supersonic flow
if M < 1
    CLa = CLa_sub;
    CMa = -CLa_sub*SM;
else
    CLa = CLa_super;
    CMa = -CLa_super*SM;
end
CM0 = 0; % Placeholder, CM0 comes from airfoil
%--------------------------------------------------------------------------
K_b   = 1;        % Roskam VI Figure 8.52 (placeholder)
cl_dfrac = 0.625; % Roskam VI Figure 8.15 (placeholder)
cl_dth = 4.2;     % Roskam VI Figure 8.14(placeholder)
k_prime = 0.7;    % Roskam VI Figure 8.13(placeholder)
a_dfrac = 1.3;    % Roskam VI Figure 8.53(placeholder)
a_de  = K_b*cl_dfrac*cl_dth*(k_prime/CLa)*a_dfrac;
Cm_ih = -CLa;
CM_de = a_de*Cm_ih; % elevator control power -------------------> need to add this
 %CM_da = 0; % effect of ailerons on pitching moment ----> need to add this
%{
fprintf('\n\n ============================================================= \n');
fprintf('\n %s Longitudinal Stability Derivatives:', flight_phase);
fprintf('\n\n CMa = %g [1/rad]', CMa);
fprintf('\n CMa = %g [1/deg]'  , CMa*pi/180);
fprintf('\n\n CL_alpha = %g [1/rad]', CLa);
fprintf('\n CL_alpha = %g [1/deg]'  , CLa*pi/180);
fprintf('\n\n elevator control power:');
fprintf('\n CM_de = %g [1/rad]', CM_de);
fprintf('\n CM_de = %g [1/deg]', CM_de*pi/180);
fprintf('\n\n ============================================================= \n');
 %}
%--------------------------------------------------------------------------
de_plot = -5:5:5;    % [deg] elevator deflection
alpha_plot = -10:10; % [deg] side slip

figure_name = sprintf('Longitudinal Stability, for %s', flight_phase);
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
for i = 1:length(de_plot)
    CM_plot(:,i) = CM0 + CMa.*alpha_plot.*pi./180 + CM_de*de_plot(i)*pi/180;
end
plot(alpha_plot, CM_plot(:,1),  'k', 'LineWidth',3);
plot(alpha_plot, CM_plot(:,2),  'b', 'LineWidth',3);
plot(alpha_plot, CM_plot(:,3),  'c', 'LineWidth',3);
hold off
xlabel('\alpha (deg)'  ,'FontSize',18);
ylabel('CM','FontSize',18);
title_string = sprintf('%s Pitching Moment vs Angle of Attack', flight_phase);
title(title_string,'FontSize',18);
legend('\delta_e = -5 [deg]','\delta_e = 0 [deg]','\delta_e = 5 [deg]');
grid on
fig_save('Figures', figure_name)
%% ========================================================================
% Lateral stability:
%--------------------------------------------------------------------------
% wing contribution:
Cl_b_basic = -1/AR;  % completely fictional equation ---> need to find this ---> use DATCOM handbook?? (via Nicolai's suggestion)
Cl_b_Lambda = -0.05; % due to sweep? ---> what equation is this?? fig. 21.10...?
Cl_b_Gamma = 0;      % no dihedral
Cl_b_wing = Cl_b_basic + Cl_b_Lambda + Cl_b_Gamma;

% vertical tail contribution:
z_v = 1; % [m] vertical distance between cg location and vertical tail ac (place holder)
var_2115 = 0.724 + (3.06*(S_VT/S_ref))/(1 + cosd(Lambda)) + 0.009*AR; % Nicolai eq. 21.15
Cl_b_VT = -CL_a_VT*var_2115*(S_VT/S_ref)*(z_v/b);

Cl_beta = Cl_b_wing + Cl_b_VT; % lateral stability derivative

beta = sqrt((1-M^2));
CL_am = 0.2;
k   = ((CL_am)*beta)/(2*pi);
cl_dfrac = 0.625; % Roskam VI Figure 8.15
cl_dth = 4.2; % Roskam VI Figure 8.14
a_da = (cl_dfrac*cl_dth)/(0.2);
Cprimel_d = (k/beta)*(0.3);
Cl_d = a_da*Cprimel_d;
Cl_da = 2*Cl_d; % aileron control power ---> use Roskam
%{
fprintf('\n\n ============================================================= \n');
fprintf('\n %s Lateral Stability Derivatives:', flight_phase);
fprintf('\n\n Cl_beta = %g [1/rad]', Cl_beta);
fprintf('\n Cl_beta = %g [1/deg]'  , Cl_beta*pi/180);
fprintf('\n\n aileron control power:');
fprintf('\n Cl_da = %g [1/rad]', Cl_da);
fprintf('\n Cl_da = %g [1/deg]', Cl_da*pi/180);
fprintf('\n\n ============================================================= \n');
%}
%--------------------------------------------------------------------------
da_plot = -5:5:5;   % [deg] rudder deflection
beta_plot = -10:10; % [deg] side slip

figure_name = sprintf('Lateral Stability, for %s', flight_phase);
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
for n = 1:length(da_plot)
    Cl_plot(:,n) = Cl_beta.*beta_plot.*pi./180 + Cl_da*da_plot(n)*pi/180;
end
plot(beta_plot, Cl_plot(:,1),  'k', 'LineWidth',3);
plot(beta_plot, Cl_plot(:,2),  'b', 'LineWidth',3);
plot(beta_plot, Cl_plot(:,3),  'c', 'LineWidth',3);
hold off
xlabel('\beta (deg)'  ,'FontSize',18);
ylabel('C_l','FontSize',18);
title_string = sprintf('%s Rolling Moment vs Sideslip', flight_phase);
title(title_string,'FontSize',18);
legend('\delta_a = -5 [deg]','\delta_a = 0 [deg]','\delta_a = 5 [deg]');
grid on
fig_save('Figures', figure_name)
%% ========================================================================
% Directional stability:
%--------------------------------------------------------------------------
% wing contribution:
x = SM*c_bar; % [m] distance between cg location and wing ac
Cn_b_wing = CL^2*(1/(4*pi*AR) - tand(Lambda)/(pi*AR*(AR + 4*cosd(Lambda)))*(cosd(Lambda) - AR/2 - AR^2/(8*cosd(Lambda)) + 6*x*sind(Lambda)/(c_bar*AR))); % Nicolai eq. 21.22

% vertical tail contribution:
Cn_b_VT = C_VT*CL_a_VT*var_2115; % Nicolai eq. 21.21

% directional stability derivative:
Cn_beta = Cn_b_wing + Cn_b_VT; % Nicolai eq. 21.20

tau = 0.6; % placeholder (Nicolai Figure 21.14) ---> this fig could be digitized and the data interpolated
Cn_dr = 0.9*CL_a_VT*C_VT*tau; % rudder control power (Nicolai eq. 21.26)
%{
fprintf('\n\n ============================================================= \n');
fprintf('\n %s Directional Stability Derivatives:', flight_phase);
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
%}
%--------------------------------------------------------------------------
dr_plot = -5:5:5;   % [deg] rudder deflection
beta_plot = -10:10; % [deg] side slip

figure_name = sprintf('Directional Stability, for %s', flight_phase);
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
for m = 1:length(dr_plot)
    Cn_plot(:,m) = Cn_beta.*beta_plot.*pi./180 + Cn_dr*dr_plot(m)*pi/180;
end
plot(beta_plot, Cn_plot(:,1),  'k', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,2),  'b', 'LineWidth',3);
plot(beta_plot, Cn_plot(:,3),  'c', 'LineWidth',3);
hold off
xlabel('\beta (deg)'  ,'FontSize',18);
ylabel('C_n','FontSize',18);
title_string = sprintf('%s Yawing Moment vs Sideslip', flight_phase);
title(title_string,'FontSize',18);
legend('\delta_r = -5 [deg]','\delta_r = 0 [deg]','\delta_r = 5 [deg]'),%'\delta_r = 6 [deg]','\delta_r = 8 [deg]','\delta_r = 10 [deg]','Location', 'best');
grid on
fig_save('Figures', figure_name)
%% ========================================================================
end
