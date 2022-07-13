%%%LINEARIZED AND NON-LINEARIZED COMPRESSION SYSTEM
%%%DEVELOPED BY : MICHAEL NATO ANDRETY
%%%STUDENT ID : 18/431068/TK/47661

%%Surge Modelling
x0=[0;0];
t_a=linspace(0,5,2500);
[t_a,x]=ode45(@compressor,t_a,x0);
[t_a,xl]=ode45(@Linearized,t_a,x0);
figure(9)
title('Pressure Rise (\psi_p) and Mass flow Rate (\phi_c)')
hold on
subplot(4,1,1)
plot(t_a,x(:,1),'linewidth',1.5)
title('Pressure Rise(\psi_p) and Mass flow Rate(\phi_c)')
xlabel('time (s)')
ylabel('\psi_p')
subplot(4,1,2)
plot(t_a,x(:,2),'linewidth',1.5)
xlabel('time (s)')
ylabel('\phi_c')
subplot(4,1,3)
plot(t_a,xl(:,1),'linewidth',1.5)
title('Linearized Pressure Rise(\xi_2) and Mass flow Rate(\xi_1)')
xlabel('time (s)')
ylabel('\xi_2')
subplot(4,1,4)
plot(t_a,xl(:,2),'linewidth',1.5)
xlabel('time (s)')
ylabel('\xi_1')

figure(10)
title('Characteristic Curve Compressor')
hold on
plot(x(:,2),x(:,1),'linewidth',1.5)
xlabel('\phi_c')
ylabel('\psi_p')

fprintf('\n')
fprintf('Non-Linear Dynamic Equation : ')
fprintf('\n')
fprintf('y1=wH/Bg*(x2-cth*uth*x1^0.5)')
fprintf('\n')
fprintf('y2=Bg*wH*(psic0+H*(1+3/2*(x2/W-1)-1/2*(x2/W-1)^3)+Po1*2/po1*kcl*input-x1)')
fprintf('\n')
fprintf('\n')
fprintf('Linearized Dynamic Equation : ')
fprintf('\n')
fprintf('y1=wH/Bg*(x2-cth*uth*(x1+psi_eq)^0.5+phi_eq);')
fprintf('\n')
fprintf('y2=Bg*wH*(psic0+H*(1+3/2*((x2+phi_eq)/W-1)-1/2*((x2+phi_eq)/W-1)^3)+kcl*input-x1-psi_eq);')
fprintf('\n')
fprintf('\n')

%%%Surge Dynamic Equation
function dxdt=compressor(t,x)
%Define Variables
mc=0.108082;  %Predicted Steady-State Mass Flow

%%Characteristic Curve Fitting 
psic0=1.5131; %Predicted Parameters psi_0             
H=0.0005;     %Predicted Parameters H
W=0.019;      %Predicted Parameters 0
PSI_ss=psic0+H*(1+3/2*(mc/W-1)-1/2*(mc/W-1)^3);

%%Compressor Parameters
Vp=0.127^2*pi*11+0.0762^2*pi*11;   %Plenum Volume
b2=324*10^-3;                      %Impeller Blade Height
Lc=30;                     %Compressor Duct Length
Ac=pi*(b2)^2;              %Compressor Cros-section area
U=249.8;                   %Impeller Tip Speed
Zc=0.9150;                 %Compresibility
R=1.47;                    %Universal Gas Constant
T=298.15;                  %Suction Temperature
y=1.2;                     %Adiabatic Constant
MW=19.66;                  %Specific Heat Ratio of Air
ao1=300;                   %Speed of Sound in Gas
wH=ao1*(Ac/Vp/Lc)^0.5;     %Helmholtz Frequency
Bg=U/(2*wH*Lc);            %Greitzer Stability Parameter
Po1=3100431.63;            %Suction Pressure
P1=4415728.38;             %Discharge Pressure
po1=28.7017;               %Suction Gas Density
psi=(P1-Po1)/(0.5*po1*U^2);%Pressure Steady-State

%%Clearance Effects
psi_ss=0.5*po1/Po1*U^2*PSI_ss+1; %Steady-State Non-Dimensional Pressure Rise
cln=1.5*10^-3;             %Nominal Clearance
cth=0.28;                  %Throttle Constant
uth=0.01;                   %Suction Massflow Representation
input=1.5*10^-6;           %Input Clearance 
k0=0.25/(1+0.25*cln/b2);   %Clearance Efficiency
%Clearance Gain
kcl=-y/(y-1)*k0/b2*psi_ss^(1/y)*(1-psi_ss^((y-1)/y));

%%Dynamic Equation
x1=x(1);
x2=x(2);
y1=wH/Bg*(x2-cth*uth*x1^0.5);
y2=Bg*wH*(psic0+H*(1+3/2*(x2/W-1)-1/2*(x2/W-1)^3)+Po1*2/po1/U^2*kcl*input-x1);
dxdt=[y1;y2];
end

%%%Surge Dynamic Equation (Linearized)
function dxdt=Linearized(t,x)
%%Surge Modelling for NG-Compressor at PT. Kaltim Methanol Industri 
%Define Variables
mc=0.108082;  %Predicted Steady-State Mass Flow

%%Characteristic Curve Fitting 
psic0=1.5131; %Predicted Parameters psi_0             
H=0.0005;     %Predicted Parameters H
W=0.019;      %Predicted Parameters 0
PSI_ss=psic0+H*(1+3/2*(mc/W-1)-1/2*(mc/W-1)^3);

%%Compressor Parameters
Vp=0.127^2*pi*11+0.0762^2*pi*11;   %Plenum Volume
b2=324*10^-3;                      %Impeller Blade Height
Lc=30;                     %Compressor Duct Length
Ac=pi*(b2)^2;              %Compressor Cros-section area
U=249.8;                   %Impeller Tip Speed
Zc=0.9150;                 %Compresibility
R=1.47;                    %Universal Gas Constant
T=298.15;                  %Suction Temperature
y=1.2;                     %Adiabatic Constant
MW=19.66;                  %Specific Heat Ratio of Air
ao1=300;                   %Speed of Sound in Gas
wH=ao1*(Ac/Vp/Lc)^0.5;     %Helmholtz Frequency
Bg=U/(2*wH*Lc);            %Greitzer Stability Parameter
Po1=3100431.63;            %Suction Pressure
P1=4415728.38;             %Discharge Pressure
po1=28.7017;               %Suction Gas Density
psi=(P1-Po1)/(0.5*po1*U^2);%Pressure Steady-State

%%Clearance Effects
psi_ss=0.5*po1/Po1*U^2*PSI_ss+1; %Steady-State Non-Dimensional Pressure Rise
cln=1.5*10^-3;             %Nominal Clearance
cth=0.28;                  %Throttle Constant
uth=0.01;                  %Suction Massflow Representation
input=1.5*10^-6;           %Input Clearance 
k0=0.25/(1+0.25*cln/b2);   %Clearance Efficiency
%Clearance Gain
kcl=-y/(y-1)*k0/b2*psi_ss^(1/y)*(1-psi_ss^((y-1)/y));

%%Finding Equilibrium Points of Compressor 
%Equilibrium Point -4.75e-06/1.303e-07 s^3 + 2.708e-07/1.303e-07 s^2 +1.972e-07/1.303e-07
phieq=roots([-4.75e-06/1.303e-07 (2.708e-07/1.303e-07-0.5*po1*U^2/Po1/cth^2/uth^2) 0 ...
    1.972e-07/1.303e-07-1]);
phi_eq=phieq(3);
psi_eq=PSI_ss;

%%Dynamic Equation
x1=x(1);
x2=x(2);
y1=wH/Bg*(x2-cth*uth*(x1+psi_eq)^0.5+phi_eq);
y2=Bg*wH*(psic0+H*(1+3/2*((x2+phi_eq)/W-1)-1/2*((x2+phi_eq)/W-1)^3)+kcl*input-x1-psi_eq);
dxdt=[y1;y2];
end