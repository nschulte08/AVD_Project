function [S,F_mdot,Pt16_Pt6] = CycleAnalysis(alpha,PRf,PRcH,Alt,M0)
%{
=====================================================================
Parametric and Performance Analysis for Low By-pass Mixed Flow Turbofan.
---------------------------------------------------------------------
Written by: Rustin Faris
Written from: "Aircraft Engine Design" by Jack D. Mattingly, Appendix H,I
Last Modified:
---------------------------------------------------------------------

=====================================================================
Engine stages:
0 = Far upstream or free stream
1 = Inlet/diffuser entry
2 = Inlet/diffuser exit fan entry
13 = Fan exit
2.5 = Low pressure compressor exit
3 = High pressure compressor exit
3.1 = Burner entry
4 = Burner exit,Nozzle vanes entry,Modedled coolant mixer 1 entry,
    High-pressure turbine entry for PRtH definition
4.1 = Nozzle vanes exit, Coolant mixer 1 exit, High-pressure turbine entry
    for PRtH definition
4.4 = High-pressure turbine exit
4.5 = Coolant mixer 2 exit,Low-pressure turbine entry
5 = Low-pressure turbine exit
6 = Core stream mixer entry
16 = Fan bypass stream mixer entry
6A = Mixer exit, Afterburner entry
7 = Afterburner exit, Exhaust nozzle entry
8 = Exhaust nozzle throat
9 = Exhaust nozzle exit
%}
%% Parametric Cycle Analysis
%-------------------------------------------------------------------%
%---------------------------CONSTANTS-------------------------------%
%-------------------------------------------------------------------%

gc = 32.174; %lbm-ft/lbf-s2 %Newtons gravitation constant 

%-------------------------------------------------------------------%
%----------------------------INPUTS---------------------------------%
%-------------------------------------------------------------------%
%---------Atmospheric Properties-------%
[~,~,T0,~,~,~,~,~,~,~] = ATMO(35000, 'E');

%---------Flight Conditions------------%
%M0 = 1.6;
% T0 = 394.10;%R
% P0 = 3.467;%psia

%---------System Parameters------------%
Beta = 0.01; %Bleed Air Fraction
Ctol = 0; %Power takeoff shaft power coefficient for low-pressure spool
Ctoh = .0150; %Power takeoff shaft power coefficient for high-pressure spool

%---------Design Limitations-----------%
hPR = 18400; %Heating Value of Fuel BTU/lbm
eps1 = .05; %Cooling Air Mass Flow #1
eps2 = .05; %Cooling Air Mass Flow #2

%---------Pressure Ratios--------------%
PRb = .950; %Burner
PRdmax = .960; %Inlet
PRMmax = .970; %Mixer due to only wall friction
PRn = .970; %Nozzle
PRAB = .950;%Afterburner

%---------Polytropic Efficiencies------%
ef = .890; %Fan
ecL = .890;%Low pressure compressor
ecH = .900;%High pressure compressor
etH = .890;%High pressure turbine
etL = .9;%%Low pressure turbine

%---------Adiabatic Efficiencies-------%
etab = .999;%Burner
etaAB = .990;%Afterburner
etamL = 1;%Mech low efficiency spool
etamH = .990;%Mech high efficiency spool
etamPL = 1;%Mech power takeoff from low pressure spool
etamPH = .990;%Mech power takeoff from high pressure spool

%---------------Design Choices---------------%
%PRf = 3.8;%Fan pressure ratio
PRcL = PRf;%Pressure ratio low pressure compressor 
%PRcH = 16;%Pressure ratio high pressure compressor
%alpha = .4;%Bypass Ratio
Tt4 = 3200;% [R] Max temperature of the high pressure turbine entry
Tt7 = 3600;% [R] Max temperature of the nozzle entry
M6 = .4;%Mach number of core stream mixer entry
P0_P9 = 1;%pressure ratio free stream to nozzle exit

%-------------------------------------------------------------------%
%---------------------------ANALYSIS--------------------------------%
%-------------------------------------------------------------------%

%--------------------Free stream-----------------------%
%T0 = 394.1;%%%%%FOR TESTING ONLY%%%%%REMOVE%%%%%
f = 0; %fuel to air ratio
[h0,Pr0,~,~,R0,Gamma_air0,a0] = FAIR1(f,T0);%h = [BTU/lbm],
V0 = M0*a0;%Inlet streamtube velocity[ft/s]
ht0 = h0+((V0.^2)/(2*778*gc));%Total enthalpy
[~,Prt0,~,~,~,~,~] = FAIR2(f,ht0);%Inlet thermal properties
taur = ht0/h0; %free stream enthalpy recovery ratio
PRr = Prt0/Pr0;%free stream pressure recovery

%---------------------Diffuser-------------------------%

%Ram recovery efficiency
if M0<=1
    etaRspec = 1;
end
if M0>1
    etaRspec = 1-0.075*(M0-1).^1.35;
end

PRd = PRdmax*etaRspec;%Pressure ratio of the diffuser
ht2 = ht0;%Adiabatic compression in diffuster
Prt2 = Prt0;%No loss in total pressure in diffuser

%----------------------FAN------------------------------%

Prt13 = Prt2*PRf.^(1/ef);%Total pressure at fan exit using defined pressure ratio for fan and defined fan polytropic efficiency
[Tt13,ht13,~,~,~,~,~] = FAIR3(f,Prt13);%Total temperature and total enthalpy based on total pressure at fan
tauf = ht13/ht2;%Enthalpy ratio accross the fan
Prt13i = Prt2*PRf;%Ideal total pressure of fan exit.
[~,ht13i,~,~,~,~,~] = FAIR3(f,Prt13i);%Ideal total enthalpy based on ideal total pressure at fan exit
etaf = (ht13i-ht2)/(ht13-ht2);%Adiabatic efficiency of the fan

%-----------------Low pressure compressor----------------%

Prt25 = Prt2*PRcL.^(1/ecL);%Total pressure at low pressure compressor exit with defined LPC pressure ratio and LPC's defined polytropic efficiency
[~,ht25,~,~,~,~,~] = FAIR3(f,Prt25);%Total enthalpy at LPC exit based on pressure at LPC exit 
taucL = ht25/ht2;%Enthalphy ratio of LPC
Prt25i = Prt2*PRcL;%Ideal total pressure at LPC exit
[~,ht25i,~,~,~,~,~] = FAIR3(f,Prt25i);%Ideal total enthalpy at LPC exit based on pressure at LPC exit
etacL = (ht25i-ht2)/(ht25-ht2);%Adiabatic efficiency of the LPC

%-----------------High pressure compressor---------------%

Prt3 = Prt25*PRcH.^(1/ecH);%Total pressure at HPC exit using polytropic efficiency
[~,ht3,~,~,~,~,~] = FAIR3(f,Prt3);%Total enthalpy at HPC exit with total pressure
taucH = ht3/ht25;%Enthalpy ratio across HPC
Prt3i = Prt25*PRcH;%Ideal total pressure at HPC exit
[~,ht3i,~,~,~,~,~] = FAIR3(f,Prt3i);%Ideal total temperature and total enthalpy at HPC exit
etacH = (ht3i-ht25)/(ht3-ht25);%Adiabatic efficiency HPC

%-----------------------Burner----------------------------%

f4i = 1;%fuel to air ratio guess for burner exit
Gate = 1;
while Gate==1
[ht4,~,~,~,~,~,~] = FAIR1(f4i,Tt4);%Ideal total enthalpy at burner exit based on fuel/air at burner exit and ideal total temperature
f = (ht4-ht3)/(etab*hPR-ht4);%defining fuel/air ratio (pg. 376, eq.6.36 in EOP)
if abs(f-f4i)>0.0001
     f4i=f;
else
    Gate = 0;
end
end

taulam = ht4/h0;%Total to static enthalpy ratio from inlet to burner exit

%------------------------Coolant mixxer 1--------------------------%
%%%%%%%%%%%%%TAUM1 and TAUTH ARE PROBLEMS MAYBE%%%%%%%%%%%%%%%%%%%%
taum1 = ((1-Beta-eps1-eps2)*(1+f)+eps1*taur*taucL*taucH/taulam)/((1-Beta-eps1-eps2)*(1+f)+eps1);%Enthalpy ratio across coolant mixer 1. (EQ 4.20a pg.111 'Aircraft Engine Design')
tautH = 1-((taur*taucL*(taucH-1)+(1+alpha)*(Ctoh/etamPH))/(etamH*taulam*((1-Beta-eps1-eps2)*(1+f)+eps1*taur*taucL*taucH/taulam)));%Enthalpy ratio across HPT (EQ 4.21a pg.112 'Aircraft Engine Design')
ht41 = ht4*taum1;%Total enthalpy at mixxer 1 exit
f41 = f/(1+f+eps1/(1-Beta-eps1-eps2));%Fuel/air at mixxer 1 exit. (EQ 4.8i pg.105 'Aircraft Engine Design')
[~,Prt41,~,~,~,~,~] = FAIR2(f41,ht41);%Total pressure at mixxer 1 exit

%--------------------High pressure turbine------------------%

ht44 = ht41*tautH;%Total enthalpy at HPT exit
[~,Prt44,~,~,~,~,~] = FAIR2(f41,ht44);%Total pressure at HPT exit
PRtH = (Prt44/Prt41).^(1/etH);%Pressure ratio at HPT exit
Prt44i=PRtH*Prt41;%Ideal total pressure at HPT exit
[~,ht44i,~,~,~,~,~] = FAIR3(f41,Prt44i);%Ideal total enthalpy at HPT exit
etatH = (ht41-ht44)/(ht41-ht44i);%Adiabatic efficiency of HPT

%----------------------Coolant mixxer 2--------------------%

taum2 = ((1-Beta-eps1-eps2)*(1+f41)+eps1+eps2*(taur*taucL*taucH/(taulam*taum1*tautH)))/((1-Beta-eps1-eps2)*(1+f41)+eps1+eps2);%Enthalpy ratio across coolant mixxer 2. (EQ 4.20b pg.112 'Aircraft Engine Design')
ht45 = ht44*taum2;%Total enthalpy at CM2 exit
f45 = f/(1+f41+(eps1+eps2)/(1-Beta-eps1-eps2));%Fuel/air at CM2 exit.(EQ 4.8j pg.105 'Aircraft Engine Design')
[~,Prt45,~,~,~,~,~] = FAIR2(f45,ht45);%Total pressure at CM2 exit

%---------------------Low pressure turbine-----------------%

tautL = 1 - (taur*((taucL-1)+alpha*(tauf-1))+(1+alpha)*Ctol/etamPL)/(etamL*taulam*tautH*((1-Beta-eps1-eps2)*(1+f)+(eps1+(eps2/tautH))*(taur*taucL*taucH/taulam)));%Enthalpy ratio across LPT. (EQ 4.22a pg.112 'Aircraft Engine Design')
ht5 = ht45*tautL;%Total enthalpy at LPT exit
[Tt5,Prt5,~,~,~,~,~] = FAIR2(f45,ht5);%Total temperature and total pressure at LPT exit
PRtL = (Prt5/Prt45).^(1/etL);%Pressure ratio across LPT
Prt5i = PRtH*Prt45;%Ideal total pressure at LPT exit
[~,ht5i,~,~,~,~,~] = FAIR3(f45,Prt5i);%Ideal total enthalpy at LPT exit
etatL = (ht45-ht5)/(ht45-ht5i);%Adiabatic efficiency of LPT

ht6 = ht5;
Tt6 = Tt5;
f6 = f45;
ht16 = ht13;
Tt16 = Tt13;
Prt16 = Prt13;
f16 = 0;

%---------------------Core stream air mixxer------------------%

alphaM = alpha/((1-Beta-eps1-eps2)*(1+f)+eps1+eps2);%Mixxer by-pass ratio (pg. 113 'Aircraft Engine Design')
f6A = f6/(1+alphaM);%Fuel/air at mixxer exit.(EQ 4.8k pg.105 'Aircraft Engine Design')
ht6A = (ht6+alphaM*ht16)/(1+alphaM);%Total enthalpy at mixxer exit
tauM = ht6A/ht6;%Enthalpy ratio at mixxer exit
Pt16_Pt6 = PRf/(PRcL*PRcH*PRb*PRtH*PRtL);%
[TSTR6,TSPR6,MFP6] = RGCOMPR1(Tt6,f45,M6);%Total to static temperature and pressure ratios and mass flow parameter at mixxer exit
T6 = Tt6/TSTR6;%Static temperature at mixxer exit

%%%%%%%%%%%%%%%ADDED(NOT IN MATTINGLY)%%%%%%%%%%%%%%%%%%%%
[~,~,~,~,R6,Gamma_air6,~] = FAIR1(f,T6);
[Tt6A,~,~,~,~,~,~] = FAIR2(f6A,ht6A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Fan by-pass air mixxer-------------------%

TSPR16 = TSPR6*Pt16_Pt6;%Total to static pressure ratio at mixxer exit 
Pr16 = Prt16/TSPR16;%Pressure at mixxer exit
[Tt16,h16,~,~,~,Gamma_air16,a16] = FAIR3(f16,Pr16);%Total temperature,total enthalpy,ratio of specific heats and speed of sound at mixxer exit
V16 = sqrt(2*gc*778*(ht16-h16));%Velocity at mixxer exit
M16 = V16/a16;%Mach number at mixxer exit
[~,~,MFP16] = RGCOMPR1(Tt16,f16,M16);%mass flow parameter at mixxer exit
A16_A6 = alphaM*sqrt(Tt16/Tt6)*(Pt16_Pt6^-1)*(MFP6/MFP16);%Area ratio of by-pass duct to mixxer exits
A6_A6A = 1/(1+A16_A6);%Area ratio of core stream mixxer to fan by-pass air mixxer exits
Constant = sqrt(R6*T6/Gamma_air6)*((1+Gamma_air6*M6.^2)+A16_A6*(1+Gamma_air16*M16.^2))/(M6*(1+alphaM));

M6Ai = .1;%initial guess of Mach number at mixxer exit

Gate = 1;
while Gate==1
[TSTR6A,~,MFP6A] = RGCOMPR1(Tt6A,f6A,M6Ai);
T6A = Tt6A/TSTR6A;
[~,~,~,~,R6A,Gamma_air6A,~] = FAIR1(f6A,T6A);
M6A = sqrt(R6A*T6A/Gamma_air6A)*((1+Gamma_air6A*M6Ai.^2)/Constant);
if abs(M6A-M6Ai)>.001
    M6Ai = M6A;
else
    Gate = 0;
end
end

PRMideal = (1+alphaM)*sqrt(tauM)*A6_A6A*(MFP6/MFP6A);%Ideal pressure ratio of the mixxer
PRM = PRMmax*PRMideal;%Pressure ratio of mixxer is pressure ratio ideal times pressure ratio due to only wall friction


f7i = .5;%initial guess of fuel/air for afterburner

Gate = 1;
while Gate==1
    [ht7,Prt7,~,~,~,~,~] = FAIR1(f7i,Tt7);
    taulamAB = ht7/h0;
    fAB = (1+f*((1-Beta-eps1-eps2)/(1+alpha-Beta)))*((taulamAB-taulam*taum1*tautH*taum2*tautL*tauM)/(hPR*etaAB/(h0-taulamAB)));
    f7 = f6A+fAB;
    if abs(f7-f7i)>.0001
        f7i=f7;
    else Gate = 0;
    end
end

%--------------------Nozzle----------------------------%

f0 = f7;
Tt9 = Tt7;%Adiabatic nozzle
ht9 = ht7;%Adiabatic nozzle
Prt9 = Prt7;%No pressure losses in nozzle

%--------------------Nozzle exit-----------------------%

TSPR9 = (P0_P9)*PRr*PRd*PRcL*PRcH*PRb*PRtH*PRtL*PRM*PRAB*PRn;%Total to static pressure ratio at nozzle exit
Pr9 = Prt9/TSPR9;%Pressure at nozzle exit
[T9,h9,~,~,R9,~,a9] = FAIR3(f0,Pr9);%Total temperature,total enthalpy,gas constant and speed of sound at nozzle exit
V9 = sqrt(2*778*gc*(ht9-h9));%Velocity at nozzle exit
M9 = V9/a9;%Mach number at nozzle exit
F_mdot = (a0/gc)*((1+f0-Beta/(1+alpha))*(V9/a0)-M0+(1+f0-Beta/(1+alpha))*((R9*(T9/T0)*(1-(P0_P9)))/(R0*(V9/a0)*Gamma_air0)));%Specific thrust at nozzle exit
S = (f0/(F_mdot))*60^2;%Uninstalled thrust specific fuel consumption
etaP = ((2*gc*M0/a0)*F_mdot)/((1+f0-(Beta/(1+alpha)))*(V9/a0).^2-M0.^2);%Uninstalled propulsive efficiency
etaTH = (1/(2*778*gc))*(((1+f0-(Beta/(1+alpha)))*V9.^2-V0.^2)+(Ctol+Ctoh)*h0)/(f0*hPR);%Uninstalled thermal efficiency
eta0 = etaTH*etaP;%Uninstalled total efficiency

%-------------------------------------------------------------------%
%---------------------------OUTPUTS---------------------------------%
%-------------------------------------------------------------------%
% Pt16_Pt6
% F_mdot
% S
% f0
% etaP
% etaTH
% V9_a0 = V9/a0
% TSPR9
% PRtH 
% PRtL 
% PRM
% tauf
% taucL
% taucH
% tautH
% tautL
% taulam
% taulamAB
% f
% fAB
% etaf
% etacL
% etacH
% etatH
% etatL

end
