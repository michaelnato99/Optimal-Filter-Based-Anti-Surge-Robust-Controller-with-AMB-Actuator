%%3D-PLOT : RADIAL MAGNETIC BEARING AND ROTOR SYSTEM

figure(3)
subplot(1,2,1)
plot(out.xcnc.data,out.ycnc.data)
title('Unconntrolled Rotor Position at Compressor Side')
xlabel('x_A (m)')
ylabel('y_A (m)')
axis(0.8*[-10^-5 10^-5 -10^-5 10^-5])

subplot(1,2,2)
plot(out.xtnc.data,out.ytnc.data)
title('Unconntrolled Position at Turbine Side (Uncontrolled)')
xlabel('x_B (m)')
ylabel('y_B (m)')
axis(0.8*[-10^-5 10^-5 -10^-5 10^-5])

figure(4)
subplot(1,2,1)
plot(out.xc.data,out.yc.data)
title('Controlled Rotor Position at Compressor Side')
xlabel('x_A (m)')
ylabel('y_A (m)')
axis(0.8*[-10^-5 10^-5 -10^-5 10^-5])

subplot(1,2,2)
plot(out.xt.data,out.yt.data)
title('Controlled Rotor Position at Turbine Side')
xlabel('x_B (m)')
ylabel('y_B (m)')
axis(0.8*[-10^-5 10^-5 -10^-5 10^-5])