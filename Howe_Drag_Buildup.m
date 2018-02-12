%%This script performs the drag buildup for subsonic, transonic, and
%%supersonic flight regimes
clear all
%DATCOM INPUTS:
%M_sub = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95];
M_sup = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8];
M_sub = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95];
%CL_sub = [0.5876, 0.33052, 0.21153, 0.1469, 0.10793, 0.08263, 0.06529, 0.0586];
CL_sup = [0.05288, 0.04371, 0.03672, 0.03129, 0.02698, 0.0235, 0.02066, 0.0183, 0.01632];
CL_sub = [0.21153, 0.1469, 0.10793, 0.08263, 0.06529, 0.0586];

%% Subsonic Drag Buildup (0<Mcr<.95)
%Inputs
pcl = 0.8;                       %Placeholder, percentage of chord that has laminar flow
Rw = 5.5;                       %Placeholder, ratio of Swet/Sref (5.5 from Howe book)
Af = .93;                       %Placeholder for airfoil factor (.93 for high-tech)
t_c = .03;                      %Placeholder 3 percent thickness for transonic (from patent)
Tf = 1.1;                       %Type factor, 1.1 for streamlined vehicle
AR = 2.7563;                    %AS2 aspect ratio
Ne = 3;                         %Number of engines
TR = 0.43847;                   %Taper Ratio
S = 125.42;                     %Sref (in m^2 from Figure 6.1 Howe)

%Drag Calc's
tau = (Rw-2)/Rw + 1.9/Rw*(1+.526*(t_c/.25)^3);      %should be close to 1
f_lam = .005*(1 + 1.5*(TR - .6)^2);                   %Taper Ratio Function, close to .0062
CDz_sub = 0.005.*(1-(2*pcl/Rw))*tau.*(1 - 0.2.*M_sub + 0.12.*(M_sub/(Af - t_c)).^20)*Rw*Tf*S^-0.1;
CDi_sub = CL_sub.^2.*((1+.12*M_sub.^6)./(pi*AR)*(1+(.142 + f_lam*AR*(10*t_c)^.33) + 0.1*(3*Ne+1)/(4+AR)^.8));

CD_sub = CDz_sub + CDi_sub;
figure(1)
plot(M_sub, CD_sub); hold on
xlabel('Mach Number'); ylabel('C_{D}');
% figure(2)
% plot(Mcr_sub, CDz_sub); hold on
% xlabel('Mach Number'); ylabel('C_{D0}');

%% Supersonic Drag Buildup
%Inputs
l = 51.8;                       %Length in meters
d = 2.635;                      %Max diameter in meters
S_l = S/(l^2);                  %Term in eq
fine = l/d;                     %Fineness/Slenderness Ratio
Ko = 1.25;                       %1.0 for ideal area distribution, 1.25 for SS airliner
Kf = 1.7;                       %1.7 for Ko = 1.25
Kw = 0.2;                       %From Howe for uncambered airfoil

%Drag Calc's
CDw = Ko/(fine^2)*(9.4*Kf/((fine^2)*(S_l)) + 1.2*((Kf/AR)^.5)*(S_l^.5)*(t_c/.05)).*(1 + .034.*(3 - M_sup).^3.5);
CDz_sup = .005.*(1 - 0.2.*M_sup)*tau*Rw*Tf*S^-.1 + CDw;
CDi_sup = CL_sup.^2.*(0.24/AR + Kw.*(M_sup.^2 - 1).^0.5);

CD_sup = CDz_sup + CDi_sup;
figure(1)
plot(M_sup, CD_sup); hold on
xlabel('Mach Number'); ylabel('C_{D}');
% figure(2)
% plot(Mcr_sup, CDz_sup); hold on
% xlabel('Mach Number'); ylabel('C_{D0}');
