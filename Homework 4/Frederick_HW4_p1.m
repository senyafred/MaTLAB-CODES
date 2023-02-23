%  stability regions of the Modified Euler method
% Z = x +1i *y , 
clear all;

ZR = -3:0.01:3;   %Real
ZI = -3:0.01:3;     %Imaginary

[X,Y] = meshgrid(ZR,ZI);    %Matrix
z = X+1i*Y;

%The Stabilty Region
Rho = abs(1 + z + (1/2)*z.^2);

%The command contour is used to plot the line curve
contour(X,Y,Rho,[1,1], "k" , "linewidth" ,2);

title( 'Stability region of Modified Euler')
xlabel( ' Re( h\lambda)')
ylabel( 'Im( h\lambda)')
axis square , grid on
xline(0); yline(0);