clear; clc; close all;

tic;

AR = 8;
TR = 0.3;
tc = 0.16;

x = 1;

for p = 1:length(tc)
    for m = 1:length(AR)
        for n = 1:length(TR)
            [OUTPUT(:,x)] = synth_config_trade(AR(m), TR(n), tc(p));
            x = x+1;
        end
    end
end

%{
OUTPUT = {num_pass; alt_cr_sub; M_cr_sub; alt_cr_super; M_cr_super;...
          AR_unswept; TR; tcmax; tmax; SM; ne;...
          WingLoading; ThrustLoading;...
          Sref; b_unswept; c_r; c_t;...
          MTOW/1000; W_fuel/1000; W_empty/1000; W_land/1000;...
          R_total_sub/1000; dt_total_sub/3600; R_total_super/1000; dt_total_super/3600;...
          V_stall; V_TO; BFL; V_approach; V_TD; FAR_land;...
          sweep_deg_TO; sweep_climb_super(2) ; sweep_deg_cr_sub ;sweep_deg_cr_super ; sweep_descend_super(1); sweep_deg_Land;...
          RTDE_Cost; DOC_Cost;...
          AB_sub; AB_super};
%}

toc;