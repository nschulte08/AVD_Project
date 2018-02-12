function [TSFC,F_mdot] = CycleAnalysis(alpha,PRf,PRcH,TechLevel)
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
clear
clc
close all

%-------------------------------------------------------------------%
%---------------------------CONSTANTS-------------------------------%
%-------------------------------------------------------------------%

gc = 32.174; %Newtons gravitation constant
% Gamma_air = 1.4;
% R = 0;
% a = 0; %Speed of sound



%-------------------------------------------------------------------%
%----------------------------INPUTS---------------------------------%
%-------------------------------------------------------------------%
%---------Atmospheric Properties-------%
[h0, P0, T0, rho0, s0, mu0, nu0, DELTA0, THETA0, SIGMA0] = ATMO(35000, 'E');

%---------Flight Conditions------------%
M0 = 1.6;
% T0 = 394.10;%R
% P0 = 3.467;%psia

%---------System Parameters------------%
Beta = 0; %Bleed Air Fraction
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
PRAB = 1;%Afterburner

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
etamPL = .995;%Mech power takeoff from low pressure spool
etamPH = .995;%Mech power takeoff from high pressure spool

%---------------Design Choices---------------%
PRf = 3.8;%Pressure ratio fan
PRcL = PRf;%Pressure ratio low pressure compressor
PRcH = 1;%Pressure ratio high pressure compressor
alpha = .4;%Bypass Ratio
Tt4 = 3200;%Max temperature of the high pressure turbine entry
Tt7 = 3600;%Max temperature of the nozzle entry
M6 = .8;%Mach number of core stream mixer entry
P0_P9 = 1;%pressure ratio free stream to nozzle exit

%-------------------------------------------------------------------%
%---------------------------ANALYSIS--------------------------------%
%-------------------------------------------------------------------%

%--------------------Free stream-----------------------%

f = 0; %fuel to air ratio
[h0,Pr0,phi0,cp0,R0,Gamma_air0,a0] = FAIR1(f,T0);
V0 = M0*a0;%Inlet velocity
ht0 = h0+(V0.^2)/(2*gc);%Total enthalpy
[Tt0,Prt0,phit0,cpt0,Rt0,Gamma_airt0,at0] = FAIR2(f,ht0);%Inlet thermal properties
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
[Tt13,ht13,phit13,cpt13,Rt13,Gamma_airt13,at13] = FAIR3(f,Prt13);%Total temperature and total enthalpy based on total pressure at fan
tauf = ht13/ht2;%Enthalpy ratio accross the fan
Prt13i = Prt2*PRf;%Ideal total pressure of fan exit.
[Tt13i,ht13i,phit13i,cpt13i,Rt13i,Gamma_airt13i,at13i] = FAIR3(f,Prt13i);%Ideal total enthalpy based on ideal total pressure at fan exit
etaf = (ht13i-ht2)/(ht13-ht2);%Adiabatic efficiency of the fan

%-----------------Low pressure compressor----------------%

Prt25 = Prt2*PRcL.^(1/ecL);%Total pressure at low pressure compressor exit with defined LPC pressure ratio and LPC's defined polytropic efficiency
[Tt25,ht25,phit25,cpt25,Rt25,Gamma_airt25,at25] = FAIR3(f,Prt25);%Total enthalpy at LPC exit based on pressure at LPC exit 
taucL = ht25/ht2;%Enthalphy ratio of LPC
Prt25i = Prt2*PRcL;%Ideal total pressure at LPC exit
[Tt25i,ht25i,phit25i,cpt25i,Rt25i,Gamma_airt25i,at25i] = FAIR3(f,Prt25i);%Ideal total enthalpy at LPC exit based on pressure at LPC exit
etacL = (ht25i-ht2)/(ht25-ht2);%Adiabatic efficiency of the LPC

%-----------------High pressure compressor---------------%

Prt3 = Prt25*PRcH.^(1/ecH);%Total pressure at HPC exit using polytropic efficiency
[Tt3,ht3,phit3,cpt3,Rt3,Gamma_airt3,at3] = FAIR3(f,Prt3);%Total enthalpy at HPC exit with total pressure
taucH = ht3/ht25;%Enthalpy ratio across HPC
Prt3i = Prt25*PRcH;%Ideal total pressure at HPC exit
[Tt3i,ht3i,phit3i,cpt3i,Rt3i,Gamma_airt3i,at3i] = FAIR3(f,Prt3i);%Ideal total temperature and total enthalpy at HPC exit
etacH = (ht3i-ht25)/(ht3-ht25);%Adiabatic efficiency HPC

%-----------------------Burner----------------------------%

f4i = .5;%fuel to air ratio guess for burner exit

Gate = 1;
while Gate==1
[ht4i,Prt4i,phit4i,cpt4i,Rt4i,Gamma_airt4i,at4i] = FAIR1(f4i,Tt3i);%Ideal total enthalpy at burner exit based on fuel/air at burner exit and ideal total temperature
f = (ht4i-ht3)/(etab*hPR-ht4i);%defining new fuel/air ratio
ht4 = ht4i;
if abs(f-f4i)>0.0001
     f4i=f;
else
    Gate = 0;
end
end

taulam = ht4/h0;%Total to static enthalpy ratio from inlet to burner exit

%------------------------Coolant mixxer 1--------------------------%

taum1 = (1-Beta-eps1-eps2)*(1+f)+eps1*taur*taucL*taucH/taulam/((1-Beta-eps1-eps2)*(1+f)+eps1);%Enthalpy ratio across coolant mixer 1
tautH = 1-((taur*taucL*(taucH-1)+(1+alpha)*(Ctoh/etamPH))/(etamH*taulam*((1-Beta-eps1-eps2)*(1+f)+eps1*taur*taucL*taucH/taulam)));%Enthalpy ratio across HPT
ht41 = ht4*taum1;%Total enthalpy at mixxer 1 exit
f41 = f/(1+f+eps1/(1-Beta-eps1-eps2));%Fuel/air at mixxer 1 exit
[Tt41,Prt41,phit41,cpt41,Rt41,Gamma_airt41,at41] = FAIR2(f41,ht41);%Total pressure at mixxer 1 exit

%--------------------High pressure turbine------------------%

ht44 = ht41*tautH;%Total enthalpy at HPT exit
[Tt44,Prt44,phit44,cpt44,Rt44,Gamma_airt44,at44] = FAIR2(f41,ht44);%Total pressure at HPT exit
PRtH = (Prt44/Prt41).^(1/etH);%Pressure ratio at HPT exit
Prt44i=PRtH*Prt41;%Ideal total pressure at HPT exit
[Tt44i,ht44i,phit44i,cpt44i,Rt44i,Gamma_airt44i,at44i] = FAIR3(f41,Prt44i);%Ideal total enthalpy at HPT exit
etatH = (ht41-ht44)/(ht41-ht44i);%Adiabatic efficiency of HPT

%----------------------Coolant mixxer 2--------------------%

taum2 = ((1-Beta-eps1-eps2)*(1+f41)+eps1+eps2*(taur*taucL*taucH/(taulam*taum1*tautH)))/((1-Beta-eps1-eps2)*(1+f41)+eps1+eps2);%Enthalpy ratio across coolant mixxer 2
ht45 = ht44*taum2;%Total enthalpy at CM2 exit
f45 = f/(1+f41+(eps1+eps2)/(1-Beta-eps1-eps2));%Fuel/air at CM2 exit
[Tt45,Prt45,phit45,cpt45,Rt45,Gamma_airt45,at45] = FAIR2(f45,ht45);%Total pressure at CM2 exit

%---------------------Low pressure turbine-----------------%

tautL = 1 - (taur*((taucL-1)+alpha*(tauf-1))+(1+alpha)*Ctol/etamPL)/(etamL*taulam*tautH*((1-Beta-eps1-eps2)*(1+f)+(eps1+(eps2/tautH))*(taur*taucL*taucH/taulam)));%Enthalpy ratio across LPT
ht5 = ht45*tautL;%Total enthalpy at LPT exit
[Tt5,Prt5,phit5,cpt5,Rt5,Gamma_airt5,at5] = FAIR2(f45,ht5);%Total temperature and total pressure at LPT exit
PRtL = (Prt5/Prt45).^(1/etL);%Pressure ratio across LPT
Prt5i = PRtH*Prt45;%Ideal total pressure at LPT exit
[Tt5i,ht5i,phit5i,cpt5i,Rt5i,Gamma_airt5i,at5i] = FAIR3(f45,Prt5i);%Ideal total enthalpy at LPT exit
etatL = (ht45-ht5)/(ht45-ht5i);%Adiabatic efficiency of LPT

ht6 = ht5;
Tt6 = Tt5;
f6 = f45;
ht16 = ht13;
Tt16 = Tt13;
Prt16 = Prt13;
f16 = 0;

%---------------------Core stream air mixxer------------------%

alphaM = alpha/((1-Beta-eps1-eps2)*(1+f)+eps1+eps2);%Mixxer by-pass ratio
f6A = f6/(1+alphaM);%Fuel/air at mixxer exit
ht6A = (ht6+alphaM*ht16)/(1+alphaM);%Total enthalpy at mixxer exit
tauM = ht6A/ht6;%Enthalpy ratio at mixxer exit
Pt16_Pt6 = PRf/(PRcL*PRcH*PRb*PRtH*PRtL);%
[TSTR6,TSPR6,MFP6] = RGCOMPR1(Tt6,f45,M6);%Total to static temperature and pressure ratios and mass flow parameter at mixxer exit
T6 = Tt6/TSTR6;%Static temperature at mixxer exit

%---------------------Fan by-pass air mixxer-------------------%

TSPR16 = TSPR6*Pt16_Pt6;%Total to static pressure ratio at mixxer exit 
Pr16 = Prt16/TSPR16;%Pressure at mixxer exit
[Tt16,h16,phi16,cp16,R16,Gamma_air16,a16] = FAIR3(f16,Pr16);%Total temperature,total enthalpy,ratio of specific heats and speed of sound at mixxer exit
V16 = sqrt(2*gc*(ht16-h16));%Velocity at mixxer exit
M16 = V16/a16;%Mach number at mixxer exit
[TSTR16,TSPR16,MFP16] = RGCOMPR1(Tt16,f16,M16);%Total to static temperature and pressure ratios and mass flow parameter at mixxer exit
A16_A6A = alphaM*sqrt(Tt16/Tt6)*(Pt6/Pt16)*(MFP6/MFP16);%Area ratio of by-pass duct to mixxer exits
A6_A6A = 1/(1+A16_A6A);%Area ratio of core stream mixxer to fan by-pass air mixxer exits
Constant = sqrt(R6*T6/Gamma_air6)*((1+Gamma_air6*M6.^2)+A16_A6A*(1+Gamma_air16*M16.^2))/(M6*(1+alphaM));
[Tt6A,Prt6A,phit6A,cpt6A,Rt6A,Gamma_airt6A,at6A] = FAIR2(f6A,ht6A);%Total temperature at mixxer exit

M6Ai = .5;%initial guess of Mach number at mixxer exit

Gate = 1;
while Gate==1
[TSTR6A,TSPR6A,MFP6A] = RGCOMPR1(Tt6A,f6A,M6Ai);
T6A = Tt6A/TSTR6A;
[h6A,Pr6A,phi6A,cp6A,R6A,Gamma_air6A,a6A] = FAIR1(f6A,T6A);
M6A = sqrt(R6A*T6A/Gamma_air6A)*((1+Gamma_air6A*M6Ai.^2)/Constant);
if abs(M6A-M6Ai)>.0001
    M6Ai = M6A;
else Gate = 0;
end
end

PRMideal = (1+alphaM)*sqrt(tauM)*A6_A6A*(MFP6/MFP6A);%Ideal pressure ratio of the mixxer
PRM = PRMmax*PRMideal;%Pressure ratio of mixxer is pressure ratio ideal times pressure ratio due to only wall friction


f7i = .5;%initial guess of fuel/air for afterburner

Gate = 1;
while Gate==1
    [ht7,Prt7,phit7,cpt7,Rt7,Gamma_airt7,at7] = FAIR1(f7i,Tt7);
    taulamAB = ht7/h0;
    fAB = (1+f*((1-Beta-eps1-eps2)/(1+alpha-Beta)))*((taulamAB-taulam*taum1*tautH*taum2*tautL*tauM)/(hPR*etaAB/(h0-taulamAB)));
    f7 = f6A+fAB;
    if abs(f7-f7i)>.0001
        f7i=f7;
    else Gate = 0;
    end
end

%--------------------Nozzle----------------------------%

f0 = f7;%All fuel burned
Tt9 = Tt7;%No total temperature loss in nozzle
ht9 = ht7;%No enthalpy losses in nozzle
Prt9 = Prt7;%No pressure losses in nozzle

%--------------------Nozzle exit-----------------------%

TSPR9 = (P0_P9)*PRr*PRb*PRcL*PRcH*PRb*PRtH*PRtL*PRM*PRAB*PRn;%Total to static pressure ratio at nozzle exit
Pr9 = Prt9*TSPR9;%Pressure at nozzle exit
[Tt9,h9,phi9,cp9,R9,Gamma_air9,a9] = FAIR3(f0,Pr9);%Total temperature,total enthalpy,gas constant and speed of sound at nozzle exit
V9 = sqrt(2*gc(ht9-h9));%Velocity at nozzle exit
M9 = V9/a9;%Mach number at nozzle exit
F_mdot = (a0/gc)*((1+f0-Beta/(1+alpha))*(V9/a0)-M0+(1+f0-Beta/(1+alpha))*((R9*(T9/T0)*(1-(P0_P9)))/(R0*(V9/a0)*Gamma_air0)));%Specific thrust at nozzle exit
S = f0/(F_mdot);%Uninstalled thrust specific fuel consumption
etaP = ((2*gc*M0/a0)*F_mdot)/((1+f0-(Beta/(1+alpha)))*(V9/a0).^2-M0.^2);%Uninstalled propulsive efficiency
etaTH = (1/2*gc)*(((1+f0-(Beta/(1+alpha)))*V9.^2-V0.^2)+(Ctol+Ctoh)*h0)/(f0*hPR);%Uninstalled thermal efficiency
eta0 = etaTH*etaP;%Uninstalled total efficiency

%-------------------------------------------------------------------%
%---------------------------OUTPUTS---------------------------------%
%-------------------------------------------------------------------%

F_mdot
S
f0
etaP
etaTH
V9_a0 = V9/a0
TSPR9
PRtH
PRtL
PRM
tauf
taucL
taucH
tautH
tautL
taulam
taulamAB
f
fAB
etaf
etacL
etacH
etatH
etatL

%% Performance Analysis

%-------------------------------------------------------------------%
%---------------------------CONSTANTS-------------------------------%
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
%----------------------------INPUTS---------------------------------%
%-------------------------------------------------------------------%

% %-------Performance choices----------%
% %Flight paramerters
% M0 = ;
% T0 = ;
% P0 = ;
% %Throttle setting
% Tt4 = ;
% Tt7 = ;
%Exhaust nozzle setting
%  P0_P9 = ;
%-------Design constants-------------%
%Pressure ratios
% PRdmax = ;
% PRb = ;
% PRMmax = ;
% PRABR = ;
% PRn = ;

%Component adiabatic efficiencies
%etaf = ;
%etacL = ;
%etacH = ;
%etatH = ;
%etatL = ;
%etab = ;
%etamL = ;
%etamH = ;
%etamPL = ;
%etamPH = ;

%Areas
%A4 =;
%A45 = ;
%A6 = ;
%A16 = ;
%A6A = ;
%A8_AB;

%Others
%Beta = ;
%eps1 = ;
%eps2 = ;
%hPR = ;
%Ptol = ;
%Ptoh = ;

%-------Reference conditions---------%
%Flight parameters
%Component behavior
%Others
%Engine control limits

%-------------------------------------------------------------------%
%---------------------------ANALYSIS--------------------------------%
%-------------------------------------------------------------------%




%-------------------------------------------------------------------%
%---------------------------OUTPUTS---------------------------------%
%-------------------------------------------------------------------%

