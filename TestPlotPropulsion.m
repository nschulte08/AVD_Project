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

%% Plots of functions used

A_pure = [2.502005e-1,-5.1536879e-5,6.5519486e-8,-6.7178376e-12,-1.5128259e-14,7.6215767e-18,-1.4526770e-21,1.0115540e-25];
A_vit = [7.3816638e-2,1.2258630e-3,-1.3771901e-6,9.9686793e-10,-4.2051104e-13,1.0212913e-16,-1.3335668e-20,7.2678710e-25];
href_pure = -1.7558886;%Btu/lbm
href_vit = 30.58153;%Btu/lbm
phi_ref_pure = 0.0454323;%[Btu/(lbm ? °R)]
phi_ref_vit = 0.6483398;%[Btu/(lbm ? °R)]
phi_ref = 1.578437947;%%%%NOT IN MATTINGLY%%%%%DETERMINED FROM TABLE DATA%%%%%%
gc = 32.174; %lbm-ft/lbf-s2 %Newtons gravitation constant

h_pure = @(T) href_pure + A_pure(1)*T+A_pure(2)/2*T.^2+A_pure(3)/3*T.^3+A_pure(4)/4*T.^4+A_pure(5)/5*T.^5+A_pure(6)/6*T.^6+A_pure(7)/7*T.^7+A_pure(8)/8*T.^8;
h_vit = @(T) href_vit + A_vit(1)*T+A_vit(2)/2*T.^2+A_vit(3)/3*T.^3+A_vit(4)/4*T.^4+A_vit(5)/5*T.^5+A_vit(6)/6*T.^6+A_vit(7)/7*T.^7+A_vit(8)/8*T.^8;
phi_pure = @(T) phi_ref_pure + A_pure(1)*log(T)+A_pure(2)*T+A_pure(3)/2*T.^2+A_pure(4)/3*T.^3+A_pure(5)/4*T.^4+A_pure(6)/5*T.^5+A_pure(7)/6*T.^6+A_pure(8)/7*T.^7;
phi_vit = @(T) phi_ref_vit + A_vit(1)*log(T)+A_vit(2)*T+A_vit(3)/2*T.^2+A_vit(4)/3*T.^3+A_vit(5)/4*T.^4+A_vit(6)/5*T.^5+A_vit(7)/6*T.^6+A_vit(8)/7*T.^7;
cp_pure = @(T) A_pure(1)+A_pure(2)*T+A_pure(3)*T.^2+A_pure(4)*T.^3+A_pure(5)*T.^4+A_pure(6)*T.^5+A_pure(7)*T.^6+A_pure(8)*T.^7;
cp_vit = @(T) A_vit(1)+A_vit(2)*T+A_vit(3)*T.^2+A_vit(4)*T.^3+A_vit(5)*T.^4+A_vit(6)*T.^5+A_vit(7)*T.^6+A_vit(8)*T.^7;


T = [3300:20:10000];
plot(T,h_pure(T))
