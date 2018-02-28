clear
clc
close all
%Shows that the maximum thrust and lowest fuel consumption occurs when the
%pressure ratio of the fan exit and the mixer are equal

PRfi= linspace(1,5,20);
PRc = 20;
Alpha = 2;

F_mdot = zeros(size(PRfi));
S = zeros(size(PRfi));
Pt16_Pt6 = zeros(size(PRfi));
PRf = zeros(size(PRfi));


for i = 1:length(PRfi)
    try
    [Si,F_mdoti,PT16_Pt6i] = CycleAnalysis(Alpha,PRfi(i),PRc,35000,1.6,1);
    catch continue
    end
    if isreal(Si)==1 && isreal(F_mdoti)==1
       S(i) = Si;
       F_mdot(i) = F_mdoti; 
       Pt16_Pt6(i) = PT16_Pt6i;
       PRf(i) = PRfi(i);
    end
end

S = S(S~=0)
F_mdot = F_mdot(F_mdot~=0);
Pt16_Pt6 = Pt16_Pt6(Pt16_Pt6~=0);
PRf = PRf(PRf~=0);

plot(PRf,Pt16_Pt6)
title('Pt16_Pt6 vs PRf')
xlabel('Pressure Ratio of the Fan')
ylabel('Ratio of Pt16/Pt6')
hold on
plot([3,4],[1,1]);
hold off

figure
plot(PRf,S)
title('TSFC vs PRf')
xlabel('Pressure Ratio of the Fan')
ylabel('Thrust Specific Fuel Consumption [lbm/hr]')

figure
plot(PRf,F_mdot)
title('Specific Thrust vs PRf')
xlabel('Pressure Ratio of the Fan')
ylabel('Specific Thrust [lbf/lbm]')

figure
plot(Pt16_Pt6,S)
title('TSFC vs Pt16/Pt6')
xlabel('Ratio of Pt16/Pt6')
ylabel('Thrust Specific Fuel Consumption [lbm/hr]')

figure
plot(Pt16_Pt6,F_mdot)
title('Specific Thrust vs Pt16/Pt6')
xlabel('Ratio of Pt16/Pt6')
ylabel('Specific Thrust [lbf/lbm]')