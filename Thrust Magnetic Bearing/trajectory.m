%%%PRE-DESIGN ANTI-SURGE CONTROLLER FULL SIMULATION
%%%DEVELOPED BY : MICHAEL NATO ANDRETY
%%%STUDENT ID : 18/431068/TK/47661

function trajectory()
ts=[-5 5];
figure(15)
hold on

%Phase Potrait Plotting
for x10=-0.1:0.02:0.1
    for x20=-0.1:0.02:0.1
    xin=[x10;x20];
    [tout, out]=ode45(@pp, ts, xin);
    x1out=out(:,1);
    x2out=out(:,2);
    plot(x1out,x2out);
    title('Phase Potrait : Linearized Pressure Rise(x) and Mass Flow Rate(\xi)')
    xlabel('PSI_P')
    ylabel('PHI_p')
    quiver(x1out,x2out,gradient(x1out),gradient(x2out),'r')
    grid on
    end
end
end

%%%Surge Dynamic Equation for Phase Potrait
function dxdt=pp(ts,x)
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
b2=324*10^-3;              %Impeller Blade Height
Lc=0.25;                   %Compressor Duct Length
Ac=pi*(b2)^2;              %Compressor Cros-section area
U=249.8;                   %Impeller Tip Speed
Zc=0.9150;                 %Compresibility
R=1.47;                    %Universal Gas Constant
T=298.15;                  %Suction Temperature
y=1.2;                     %Adiabatic Constant
MW=19.66;                  %Specific Heat Ratio of Air
ao1=200;                   %Speed of Sound in Gas
wH=ao1*(Ac/Vp/Lc)^0.5/4;   %Helmholtz Frequency
Bg=U/(2*wH*Lc);            %Greitzer Stability Parameter
Po1=3100431.63;            %Suction Pressure
P1=4415728.38;             %Discharge Pressure
po1=28.7017;               %Suction Gas Density
psi=(P1-Po1)/(0.5*po1*U^2);%Pressure Steady-State

%%Clearance Effects
psi_ss=0.5*po1/Po1*U^2*PSI_ss+1; %Steady-State Non-Dimensional Pressure Rise
cln=3*10^-3;               %Nominal Clearance
cth=0.28;                  %Throttle Constant
uth=0.1;                   %Suction Massflow Representation
input=10^-4;               %Input Clearance 
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