%  stability regions of the Runge - Kutta method, Simple Plot
clear all;

ZR= -3:0.01:2.5;    %Real Part

Z = ZR;

% The stability region
R4 = 1 + Z + (1/2)*Z.^2 + (1/6)*Z.^3 + (1/24)*Z.^4;

%The command contour is used to plot the line curve
plot(ZR, R4,"b" , "linewidth" ,2);hold on;
plot(ZR, yline(1),"k" , "linewidth" ,2); hold off;

title( 'Stability region of  RK4')
xlabel( ' Re( h\lambda)')
ylabel( 'Im( h\lambda)')
axis square , grid on
xline(0); yline(0);