clear; clc; close all;
tic;
%% ========================================================================
h = 16000:2000:18000;
M_cr_sub = 0.7;
M_cr_super = 1.2:0.2:2.0;

x = 1;

for n = 1:length(h)
    for m = 1:length(M_cr_super)
        [OUTPUT(:,x)] = synth_mission_trade(h(n), M_cr_sub, M_cr_super(m));
        x = x+1;
    end
end
%% ========================================================================
toc;