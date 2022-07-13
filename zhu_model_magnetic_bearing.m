%THRUST MAGNETIC BEARING PRE-DESIGN
%DEVELOPED BY : MICHAEL NATO ANDRETY
%STUDENT ID : 18/431068/TK/47661

%Define Variables
%Rotor (Rigid)
m=330;      %Rotor Mass (kg)
c=100;
k=25;
N=143;      %Number of Turns
imax=10;    %Maximal Current
ib=6.5;     %Bias Current
f0=ib*N;    %Nominal Force
fmax=imax*N;%Maximal Force

%Thrust Magnetic Bearing
max_curr_den=300/10000; %Maximal Current Density
go=3*10^-3 ;    %Nomial Air Gap        
Bsat=1.2;       %Magnetic Saturation 
miu0=4*pi*10^-7;%Permeability in Air
miur=3000;      %Iron Relative Permeability
iron_c=10^7;    %Iron Conductivity

%Dimension of Thrust Magnetic Bearing
r1=122.5*10^-3;
r2=143.8*10^-3;
r3=150.0*10^-3;
d1=7.5*10^-3;
d2=28*10^-3;
d3=7.5*10^-3;
%Cross Sectional Area
A2=pi*r3^2-pi*r2^2;
A1=pi*r1^2;
%Effective Reluctances 
R1=go/(pi*r1^2*miu0);               %Eff. Rel. on Region 1
R2=log(r2/r1)/(2*pi*d1*miur*miu0);  %Eff. Rel. on Region 2
R3=go/(pi*miu0*r3^2-r2^2);          %Eff. Rel. on Region 3
R4=d2/(pi*miur*miu0*(r3^2-r2^2));   %Eff. Rel. on Region 4
R5=log(r2/r1)/(2*pi*d3*miur*miu0);  %Eff. Rel. on Region 5
R6=d2/(pi*r1^2*miur*miu0);          %Eff. Rel. on Region 6
%Total
Ro=R1+R2+R3+R4+R5+R6;

fprintf('Effective Reluctances for Element 1~6')
[R1 R2 R3 R4 R5 R6]

%Eddy Current Approximation for each Region (1-6)
c1=(1/(4*pi))*(iron_c/(miur*miu0))^0.5;
c2=log(r2/r1)/(2*pi)*(iron_c/(miur*miu0))^0.5;
c3=(2*r3^4*log(r3/r2)-3/2*r3^4+2*r3^2*r2^2-1/2*r2^4)/(2*pi*(r3^2-r2^2)^2)*(iron_c/(miur*miu0))^0.5;
c4=d2/(2*pi*r2)*(iron_c/(miur*miu0))^0.5;
c5=log(r2/r1)/(2*pi)*(iron_c/(miur*miu0))^0.5;
c6=d2/(2*pi*r1)*(iron_c/(miur*miu0))^0.5;
%Total
c=c1+c2+c3+c4+c5+c6;
fprintf('Eddy Current Approximation for Element 1~6')
[c1 c2 c3 c4 c5 c6]

%Gain
fprintf('Current and Position Gain')
Kflux=1/miu0*(1/A1+1/A2)*f0/Ro;
%Current Ki and Position Gain Kx
Ki=Kflux*N/Ro;
Kx=Kflux^2/Ro/25000;

%Fractional Order Transfer Function of TAMB
s=fotf('s');
es=Ro/(Ro+c*s^0.5); %Eddy Current Effect
ms=1/(m*s^2+c*s+k);%Shaft
hs=es*ms;           %Open Loop FOTF
thrust_tf=Ki*(es/1-Kx*es)


figure(12)
subplot(2,1,1)
bode(hs)
title('Bode Plot : Thrust Magnetic Bearing System')
subplot(2,1,2)
bode(thrust_tf)
title('Bode Plot : Thrust Magnetic Bearing System and Rotor')

figure(14)
subplot(2,1,1)
step(hs)
title('Step Response : Thrust Magnetic Bearing System')
xlabel('time (s)')
ylabel('Force (N)')
subplot(2,1,2)
step(thrust_tf)
title('Step Response : Thrust Magnetic Bearing System and Rotor')
xlabel('time (s)')
ylabel('Rotor Position (m)')

fprintf('\n')
fprintf('KALMAN FILTERING :')
fprintf('\n')
%%Process Noise Q, Measurement Noise R
x0T=0;
P0T=10^-3;
QcovT=10^-5;
RcovT=10^-4;

%%Discrete Time Model
sys_dt=ss(0.92,0.1,1,0,0.1); %Zero Order Hold

%Discrete A,B,C,D Matrices
AolT=sys_dt.A;
BolT=sys_dt.B;
ColT=sys_dt.C;
DolT=sys_dt.D;
