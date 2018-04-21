%{
===========================================================================
Verical Tail sizing function
===========================================================================
INPUTS:
b_eff   = effective span [m]
S_ref   = referrence area [m^2]
Lambda  = sweep angle [deg]
===========================================================================
OUTPUTS:
S_VT    = total required vertical tail area [m^2]
C_VT    = vertical tail volume coefficient
y_VT    = spanwise locationof vertical tail area [m]
l_VT    = distance (in x-direction) between cg location and 1/4 chord of VT mac [m]
vert_tail_plot = data for plotting S_VT vs l_VT
===========================================================================
%}
function [S_VT, l_VT, C_VT, y_VT, vert_tail_plot] = VT_size(Lambda, S_ref, b_eff)
%% ========================================================================
% vertical tail sizing:
C_VT = 0.09;                % Nicolai table 11.8 (buisness aircraft)
lVT_SVT = C_VT*b_eff*S_ref; % [m^3] product of l_VT and S_VT (Nicolai p.286)
%--------------------------------------------------------------------------
% inital estimate:
y_VT = 0.9*(b_eff/2);     % location of VT along the left wing span = 90% of half span (arbitrary)
l_VT = y_VT*sind(Lambda); % this is a function of sweep
S_VT = lVT_SVT/l_VT;
%--------------------------------------------------------------------------
% plot data:
l_VT_plot = 1:0.1:l_VT;         % [m] distance (in x-direction) between cg location and 1/4 chord of VT mac
S_VT_plot = lVT_SVT./l_VT_plot; % [m^2] VT area

vert_tail_plot = [l_VT_plot; S_VT_plot];

end
