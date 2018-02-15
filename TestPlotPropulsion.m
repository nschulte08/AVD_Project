clear
clc
close all

%USE THIS TO CHECK THE FAIR FUNCTIONS USING PAGE 818 IN 'ELEMENTS OF PROPULSION'%

T = [3300:20:3750]';
h = linspace(878.77,1009.78,length(T));h = h';
Pr = linspace(1413.4,2434,length(T));Pr = Pr';
h_out1 = [];
Pr_out1 = [];
T_out2 = [];
Pr_out2 = [];
T_out3 = [];
h_out3 = [];

for i = 1:length(T)
    [h_out,Pr_out,~,~,~,~,~] = FAIR1(0,T(i));
    h_out1 = [h_out1;h_out];
    Pr_out1 = [Pr_out1;Pr_out];
end

for i = 1:length(T)
    [T_out,Pr_out,~,~,~,~,~] = FAIR2(0,h(i));
    T_out2 = [T_out2;T_out];
    Pr_out2 = [Pr_out2;Pr_out];
end

for i = 1:length(T)
    [T_out,h_out,~,~,~,~,~] = FAIR3(0,Pr(i));
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



