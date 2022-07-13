%%%PRE-DESIGN ANTI-SURGE CONTROLLER FULL SIMULATION
%%%Radial Magnetic Bearing Design
%%%DEVELOPED BY : MICHAEL NATO ANDRETY
%%%STUDENT ID : 18/431068/TK/47661

%%Define Variables
%Rotor (Rigid)
clear all

Fdist=2920.4;       %Distortion Force
kem=7.41828e-08;    %Electromagnetic Constant
Fmax=4205.376;      %Maximum Force
Flinmax=6055.74144; %Maximum Linear Force
speed=1151;   %Rotor Speed
gm0=0.5;            %Copper Fill Factor
Ix=7.8498;          %Rotor Inertia X-axis
Iy=Ix;              %Rotor Inertia Y-axis
Iz=0.03635;         %Rotor Inertia Z-axis
m=520;              %Rotor Mass

%Sensor Position
a=-1.707;
b=1.983;
c=-2.007;
d=2.283;

%Amplifier Design
ki=582.1076537;         %Current Amplifier Constant
ks=-67146068.109;       %Stiffness Amplifier Constant
Ki=eye(4)*ki;           %ki in Matrix Form
Ks=eye(4)*ks;           %ks in Matrix Form
Bsat=1.2;               %Saturation of Magnetic Flux Density
rfreq=183.25;           %Rotor Frequency
slewmax=4842043.464;    %Maximum Slew Rate
Vamax=2421.021732;      %Maximum Power Amplifie
Imax=10;                %Maximum Coil Current
Vmax=242.1021732;       %Maximum Coil Voltage
biasratio=0.576749512;  %Bias Ratio
Ibias=5.767495124;      %Bias Current (Offset)
Na=334;                 %Number of Turns

%Gyroscopic Matrix G
G=[0 0 Iz*speed 0;0 0 0 0 ...
    ;-Iz*speed 0 0 0 ...
    ;0 0 0 0];

%Input Matrix
B=[a b 0 0;1 1 0 0 ...
    ;0 0 a b ...
    ;0 0 1 1];

%Output Matrix
C=[c 1 0 0;d 1 0 0 ...
    ;0 0 c 1 ...
    ;0 0 d 1];

%Mass Matrix
M=[Iy 0 0 0;0 m 0 0 ...
    ;0 0 Ix 0;0 0 0 m];

%Transformation Matrix
T=[a 1 0 0;b 1 0 0 ...
    ;0 0 a 1; 0 0 b 1];
Tr=inv(T);
%Negative Bearing Stiffness Matrix
Kss=(B*Ks*B')';

Ar=[zeros(4) eye(4); zeros(4) -G/M];
Br=[zeros(4);B/M];
Cr=[C zeros(4)];
Dr=zeros(4);

%State Space Form
%Define States
states={'α','β','x_SA','x_SB','x_A','x_B','y_A','y_B'};
inputs={'i_x_A','i_x_B','i_y_A','i_y_B'};
outputs={'x_A','x_B','y_A','y_B'};
rad_amb=ss(Ar,Br,Cr,Dr,'statename',states,'inputname',inputs,'outputname',outputs);
Ad=rad_amb.A;
Bd=rad_amb.B;
Cd=rad_amb.C;
Dd=rad_amb.D;

%Transfer Function
s=tf('s');
Hsi=Ki*(feedback(rad_amb,-Kss))*Tr;
H1=Hsi(1,1);
H2=Hsi(2,2);
H3=Hsi(3,3);
H4=Hsi(4,4);

figure(1)
step(rad_amb)
title('Step Response : RMB Dynamics Model')
xlabel('time (s)')
figure(2)
bode(rad_amb)
title('Bode Plot : RMB Dynamics Model')

%%Discrete Kalman Filter
dT=0.005; %Fs=10Hz
sys_dr=c2d(Hsi,dT,'zoh');
%Process Noise Q, Measurement Noise R
x0ramb=[0;0;0;0;0;0;0;0];
P0ramb=diag(100*ones(1,8));
Qramb=diag(10^-14*ones(1,8));
Rramb=diag(10^-14*ones(1,4));

%%LQG Controller Design
%Cost Function
Crmb=ctrb(Hsi)
Ormb=obsv(Hsi)
Q=eye(8)*0.7;
R=eye(4)*0.4;
[K,Prmb,e] = lqr(Hsi,Q,R,[]);
Cc=[eye(4) zeros(4,4)] ;
Nbar= -inv(Cc*inv(Ad-Bd*K)*Bd);

%LQR Gain
KLQR=Nbar;
OL=Hsi*KLQR;