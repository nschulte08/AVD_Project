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
function [VT, VT_plot] = VT_size(Lambda, S_ref, b_eff)
%% ========================================================================
% vertical tail sizing:
C_VT = (0.02 + 0.05)/2;                % Nicolai table 11.8 (buisness aircraft)
lVT_SVT = C_VT*b_eff*S_ref; % [m^3] product of l_VT and S_VT (Nicolai p.286)
%--------------------------------------------------------------------------
% inital estimate:
y_VT = 0.9*(b_eff/2);     % location of VT along the left wing span = 90% of half span (arbitrary)
l_VT = y_VT*sind(Lambda); % this is a function of sweep
S_VT = lVT_SVT/l_VT;
%--------------------------------------------------------------------------
% Tail Geometry
TR_VT       = 0.5;
AR_VT       = 2;
b_VT        = sqrt(AR_VT*S_VT);
cr_VT       = 2*b_VT/(AR_VT*(1+TR_VT));
ct_VT       = TR_VT*cr_VT;
cbar_VT     = (2/3)*cr_VT*((1+TR_VT+(TR_VT^2))/(1+TR_VT));
Z_bar       = (b_VT/6)*((1+(2*TR_VT))/(1+TR_VT));
SweepLE_VT  = atand((cr_VT-ct_VT)/b_VT);
%--------------------------------------------------------------------------
% plot data:
l_VT_plot = 1:0.1:l_VT;         % [m] distance (in x-direction) between cg location and 1/4 chord of VT mac
S_VT_plot = lVT_SVT./l_VT_plot; % [m^2] VT area

VT_plot = [l_VT_plot; S_VT_plot];
%--------------------------------------------------------------------------
% collect results:
VT.S_VT         = S_VT;
VT.l_VT         = l_VT;
VT.C_VT         = C_VT;
VT.y_VT         = y_VT;
VT.TR           = TR_VT;
VT.AR           = AR_VT;
VT.b            = b_VT;
VT.cr           = cr_VT;
VT.ct           = ct_VT;
VT.cbar         = cbar_VT;
VT.Z_bar        = Z_bar;
VT.SweepLE_VT   = SweepLE_VT;

end
