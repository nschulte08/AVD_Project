%% Control Surface Deflections=============================================
Dela  = 0.5*(Dela_l-Dela_r);
Dele  = 
Delr  = 
Cm_ih = 0                                                        %since OFW
Ih    = 0                                                        %since OFW
%% Pitching Moment=========================================================
Cm = Cm_0+(Cm_aoa*AOA)+(Cm_dele*Dele)+(Cm_ih*Ih)
%% Rolling Moment==========================================================
Cl = Cl_0+(Cl_beta*BETA)+(Cl_dela*Dela)+(Cl_delr*Delr)
%% Yawing Moment===========================================================
Cn = Cn_0+(Cn_beta*BETA)+(Cn_dela*Dela)+(Cn_delr*Delr)