%%%%PRE-DESIGN ANTI-SURGE CONTROLLER FULL SIMULATION
%%%%DEVELOPED BY : MICHAEL NATO ANDRETY
%%%%STUDENT ID : 18/431068/TK/47661

clear all; 
close all; 
clc
%%%Surge Modelling for NG-Compressor at PT. Kaltim Methanol Industri 
%%Define Variables
mc=0.108082;  %Predicted Steady-State Mass Flow

%Characteristic Curve Fitting
psic0=1.5131; %Predicted Parameters psi_0
H=0.0005;     %Predicted Parameters H
W=0.019;      %Predicted Parameters 0
PSI_ss=psic0+H*(1+3/2*(mc/W-1)-1/2*(mc/W-1)^3);

%Compressor Parameters
Vp=0.127^2*pi*11+0.0762^2*pi*11;   %Plenum Volume
b2=324*10^-3;              %Impeller Blade Height
Lc=30;                     %Compressor Duct Length
Ac=pi*(b2)^2;              %Compressor Cros-section area
U=249.8;                   %Impeller Tip Speed
Zc=0.9150;                 %Compresibility
R=1.47;                    %Universal Gas Constant
T=298.15;                  %Suction Temperature
y=1.2;                     %Adiabatic Constant
MW=19.66;                  %Molar Mass
ao1=300;                   %Speed of Sound in Gas
wH=ao1*(Ac/Vp/Lc)^0.5;     %Helmholtz Frequency
Bg=U/(2*wH*Lc);            %Greitzer Stability Parameter
Po1=3100431.63;            %Suction Pressure
P1=4415728.38;             %Discharge Pressure
po1=28.7017;               %Suction Gas Density
psi=(P1-Po1)/(0.5*po1*U^2);%Pressure Steady-State
pu=po1;                    %Pipeline Gas Density

%Clearance Effects
%%Clearance Effects
psi_ss=0.5*po1/Po1*U^2*PSI_ss+1; %Steady-State Non-Dimensional Pressure Rise
cln=1.5*10^-3;             %Nominal Clearance
cth=0.28;                  %Throttle Constant
uth=0.01;                  %Suction Massflow Representation
input=1.5*10^-6;           %Input Clearance 
k0=0.25/(1+0.25*cln/b2);   %Clearance Efficiency
%Clearance Gain
uth_eq=0.17;               %Nominal Throttle Valve Opening
%Clearance Gain
kcl=-y/(y-1)*k0/b2*psi_ss^(1/y)*(1-psi_ss^((y-1)/y));

%%%Finding Equilibrium Points of Compressor 
%Equilibrium Point -4.75e-06/1.303e-07 s^3 + 2.708e-07/1.303e-07 s^2 +1.972e-07/1.303e-07
phieq=roots([-4.75e-06/1.303e-07 (2.708e-07/1.303e-07-0.5*po1*U^2/Po1/cth^2/uth^2) 0 ...
    1.972e-07/1.303e-07-1]);
phi_eq=phieq(3);
psi_eq=PSI_ss;

%%%State-Space Representation of Linearized Model
a1=-4.75e-06/1.303e-07;
b1=2.708e-07/1.303e-07;
a2=1.1925;
syms s               %Laplace Operator
d=2.83*10^-5;        %Line Dissipation Number
lambda1=pi*(1-0.5)/d;%Undamped Natural Frequency
Z=Zc/1000;           %Line Impedance Constant

%%Piping Dynamics Equation of Compression System
%%xn^dot = Aixn + Biu
%%y=Cxn
%Modal Aproximation
Ai=[0 (-1)^2*Z*lambda1;-((-1)^2*lambda1)/Z 8];
Bi=[0 -2*Z/d;2/(Z*d) 0];
Hn=-inv(Ai)*Bi;
G=inv(Hn)*[1 -8*Z*d;0 1];

%State Space Representation
Apipe=[0 -Z*pi/(2*d);pi/(2*d*Z) -8];
Bpipe=[0 -2*Z/d;2*Z/d 0]*G;

%%Define Complete State Space Compression System
A=[Bg*wH*(3*a1*phi_eq^2+2*b1*phi_eq) -Bg*wH 0 0;
    wH/Bg -wH/Bg*(cth*uth*a2) 0 -wH/Bg;
    0 0 (Bpipe(1,2)*Ac*uth_eq*cth)/(pu*(psi_eq^0.5)) 2*Apipe(1,2)*Ac/(pu*U);
    0 Bpipe(2,1)*pu*U/(2*Ac) (Apipe(2,1)*pu*U/(2*Ac))+(Bpipe(2,2)*uth*cth/(2*psi_eq^0.5)) Apipe(2,2)
    ];
B=[2*Bg*wH*Po1*kcl/(po1*U^2);0;0;0];
C=[1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
D=[0;0;0;0];
%Define States
states={'\xi_1','\xi_2','\xi_3','\xi_4'};
inputs={'delta'};
outputs={'\xi_1','\xi_2','\xi_3','\xi_4'};
sysc=ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

figure(1)
step(sysc)
title('Complete Linearized Compression System Model (Normal)')
xlabel('time (s)')
Co = ctrb(sysc);
Oo=obsv(sysc);
s=tf('s');
Hs=sysc.C*inv(s*eye(4)-sysc.A)*sysc.B;

%%%Hinfinity Loop-Shaping Design
dm1 = ureal('dm',2,'Percentage',5);
dn1 = ureal('dn',4,'Percentage',10);
dm2 = ureal('dm',3,'Percentage',3);
dn2 = ureal('dn',1.5,'Percentage',12);

figure(2)
subplot(2,2,1)
[fact1,Ml1,Nl1] = lncf(Hs(1));
Ml1;
Nl1;
sys1=(Nl1+dn1)/(Ml1+dm1);
sigma(sys1,'g',{.1,100000});
title('Singular value plot for Hs_1')

subplot(2,2,2)
Gd1 =tf([0 5.52683216095313 622.867528560800 17983717127.2993 1782362086895.86],...
    [1 132.359638468200 3253895914.18623 386466689536.296 10361510011795.3]);     
sigma(Gd1,{.1 100})
grid
title('Target loop shape Gd1(s).')
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')

[K1,CL1,GAM1] = loopsyn(sys1,Gd1);
GAM1
K1

subplot(2,2,3)
[fact2,Ml2,Nl2] = lncf(Hs(2));
Ml2
Nl2
sys2=(Nl2+dn2)/(Ml2+dm2)
sigma(sys2,'g',{.1,100000});
title('Singular value plot for Hs_1')

subplot(2,2,4)
Gd2 =tf([0 0 1732.63599361523 13900.4983909958 5337946481127.45],...
    [1 132.359638468200 3253895914.18623 386466689536.296 10361510011795.3]);
sigma(Gd2,{.1 100000})
grid
title('Target loop shape Gd2(s).')
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')

[K2,CL2,GAM2] = loopsyn(sys2,Gd2);
GAM2
K2

L1 = sys1*K1;              %form the compensated loop L1

figure(3)
subplot(2,1,1)
sigma(Gd1,'b',L1,'r--',{.1,100000});
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
grid on
legend('Gd1 (target loop shape)','L1 (actual loop shape)');

L2 = sys2*K2;              %form the compensated loop L2
subplot(2,1,2)
sigma(Gd2,'b',L2,'r--',{.1,100000});
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
grid on
legend('Gd2 (target loop shape)','L2 (actual loop shape)');

%Analyzing the Shaped-Loop L1, Closed-Loop T1, and Sensitivity S1
T1 = feedback(L1,eye(1));
T1.InputName = {'Clearance'};
S1 = eye(1)-T1;
%Analyzing the Shaped-Loop L2, Closed-Loop T2, and Sensitivity S2
T2 = feedback(L2,eye(1));
T2.InputName = {'Clearance'};
S2 = eye(1)-T2;

% SIGMA frequency response plots
figure(4)
subplot(2,1,1)
sigma(inv(S1),'m',T1,'g',L1,'r--',Gd1,'b',Gd1/GAM1,'b:',...
	Gd1*GAM1,'b:',{.1,100000})
grid on
legend('1/\sigma(S1) performance',...
	'\sigma(T1) robustness',...
	'\sigma(L1) open loop',...
	'\sigma(Gd1) target loop shape',...
	'\sigma(Gd1) \pm GAM1(dB)');
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
% Make lines wider and fonts larger
set(findobj(gca,'Type','line','-not','Color','b'),'LineWidth',2);
h = findobj(gca,'Type','line','-not','Color','b');
set(h,'LineWidth',2);

subplot(2,1,2)
sigma(inv(S2),'m',T2,'g',L2,'r--',Gd2,'b',Gd2/GAM2,'b:',...
	Gd2*GAM2,'b:',{.1,100000})
grid on
legend('1/\sigma(S2) performance',...
	'\sigma(T2) robustness',...
	'\sigma(L2) open loop',...
	'\sigma(Gd2) target loop shape',...
	'\sigma(Gd2) \pm GAM(dB)');
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
% Make lines wider and fonts larger
set(findobj(gca,'Type','line','-not','Color','b'),'LineWidth',2);
h = findobj(gca,'Type','line','-not','Color','b');
set(h,'LineWidth',2);

figure(5)
subplot(2,1,1)
step(T1,'b',sys1,'r--')
legend('Controlled Output \xi_1 (s)/U(s)','Real Output \xi_1 (s)/U(s)');

subplot(2,1,2)
step(T2,'b',sys2,'r--')
legend('Controlled Output \xi_2 (s)/U(s)','Real Output \xi_2 (s)/U(s)');
Mls1=Ml1.C*inv(s*eye(4)-Ml1.A)*Ml1.B;
Mls2=Ml2.C*inv(s*eye(4)-Ml2.A)*Ml2.B;

figure(6)
subplot(2,1,1)
sigma(T1,'b',Hs(1),'r--')
legend('Controlled Output \xi_1 (s)/U(s)','Real Output \xi_1 (s)/U(s)');

subplot(2,1,2)
sigma(T2,'b',Hs(2),'r--')
legend('Controlled Output \xi_2 (s)/U(s)','Real Output \xi_2 (s)/U(s)');
Mls1=Ml1.C*inv(s*eye(4)-Ml1.A)*Ml1.B;
Mls2=Ml2.C*inv(s*eye(4)-Ml2.A)*Ml2.B;

fprintf('eig(Ml1.A) :')
fprintf('\n')
lambda_Ml1=eig(Ml1.A)
fprintf('eig(Nl1.A) :')
fprintf('\n')
lambda_Nl1=eig(Nl1.A)
fprintf('eig(Ml2.A) :')
fprintf('\n')
lambda_Ml2=eig(Ml2.A)
fprintf('eig(Nl2.A) :')
fprintf('\n')
lambda_Nl2=eig(Nl2.A)

figure(7)
subplot(2,1,1)
step(sys1)
title('Step Response : Uncertainty Model of \xi_1 (s)/U(s)')
xlabel('time (s)')
subplot(2,1,2)
margin(sys1)
title('Bode Plot : Uncertainty Model of \xi_1 (s)/U(s)')

figure(8)
subplot(2,1,1)
step(sys2)
title('Step Response : Uncertainty Model of \xi_2 (s)/U(s)')
xlabel('time (s)')
subplot(2,1,2)
margin(sys1)
title('Bode Plot : Uncertainty Model of \xi_2 (s)/U(s)')

%%%Call Other Function 
surge_remodel()
fprintf('\n')
fprintf('Equilibrium Point :')
fprintf('\n')
phi_eq
psi_eq
trajectory()

fprintf('------------------------------')
fprintf('\n')
fprintf('THRUST MAGNETIC BEARING SIDE :')
fprintf('\n')
fprintf('------------------------------')
fprintf('\n')
fprintf('\n')

zhu_model_magnetic_bearing()

%%%Discrete Kalman Filter Startup
dT=0.005; %Sampling time, Fs=20Hz
sys_d=c2d(sysc,dT,'zoh'); 

figure(13)
subplot(2,1,1)
step(sys_d)
title('Step Response : Discretized Compression System')
xlabel('time (s)')
subplot(2,1,2)
pzmap(sys_d)
title('Pole & Zero Map: Discretized Compression System')

%%Process Noise Q, Measurement Noise R
x0c=[0;0;0;0];
P0=diag(0.0001^2*ones(1,4));
Qcov=diag(0.001^2*ones(1,4));
Rcov=diag(0.1^2*ones(1,4));