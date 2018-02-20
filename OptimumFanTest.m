clear
clc
close all
%Shows that the maximum thurst and lowest fuel consumption occurs when the
%pressure ratio of the fan exit and the mixer are equal

PRfi= linspace(3,4,50);
PRc = 16;
Alpha = .4;

F_mdot = [];
S = [];
Pt16_Pt6 = [];
PRf = [];


for i = 1:length(PRfi)
    [Si,F_mdoti,PT16_Pt6i] = CycleAnalysis(Alpha,PRfi(i),PRc,35000,1.6);
    if isreal(Si)==1 && isreal(F_mdoti)==1
       S = [S,Si];
       F_mdot = [F_mdot,F_mdoti];
       Pt16_Pt6 = [Pt16_Pt6,PT16_Pt6i];
       PRf = [PRf,PRfi(i)]; 
    end 
end
    
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