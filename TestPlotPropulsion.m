clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USE THIS TO CHECK THE FAIR FUNCTIONS USING PAGE 818 IN 'ELEMENTS OF PROPULSION'%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = [3300:20:4000]';
h = linspace(964.60,1202.13,length(T));h = h';%don't forget to change these if you change fuel to air or temperature range
Pr = linspace(2795,7227,length(T));Pr = Pr';
h_out1 = [];
Pr_out1 = [];
T_out2 = [];
Pr_out2 = [];
T_out3 = [];
h_out3 = [];

for i = 1:length(T)
    [h_out,Pr_out,~,~,~,~,~] = FAIR1(0.0676,T(i));%other input is fuel to air ratio
    h_out1 = [h_out1;h_out];
    Pr_out1 = [Pr_out1;Pr_out];
end

for i = 1:length(T)
    [T_out,Pr_out,~,~,~,~,~] = FAIR2(0.0676,h(i));
    T_out2 = [T_out2;T_out];
    Pr_out2 = [Pr_out2;Pr_out];
end

for i = 1:length(T)
    [T_out,h_out,~,~,~,~,~] = FAIR3(0.0676,Pr(i));
    T_out3 = [T_out3;T_out];
    h_out3 = [h_out3;h_out];
end

FAIR1 = table(T,h_out1,Pr_out1);
FAIR1.Properties.Description = 'Outputs of FAIR1 Function from input temperatures';
FAIR1.Properties.VariableUnits = {'R' 'BTU/lbm' ''};

FAIR2 = table(T_out2,h,Pr_out2);
FAIR2.Properties.Description = 'Outputs of FAIR2 Function from input enthalpys';
FAIR2.Properties.VariableUnits = {'R' 'BTU/lbm' ''};

FAIR3 = table(T_out3,h_out3,Pr);
FAIR3.Properties.Description = 'Outputs of FAIR3 Function from input reduced pressures';
FAIR3.Properties.VariableUnits = {'R' 'BTU/lbm' ''};

FAIR1
FAIR2
FAIR3



