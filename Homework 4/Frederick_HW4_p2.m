%  stability regions of the Runge - Kutta method
% Z = x +1i *y , .
clear all;

ZR = -3:0.01:3;   %Real
ZI = -3:0.01:3;     %Imaginary

[X , Y ] = meshgrid (ZR, ZI);   %Matrix
 
Z = X +1i * Y ;
 
%The Stabilty Region
R4 = abs(1 + Z + .5*Z.^2 + (1/6)*Z.^3 + (1/24)*Z.^4);

% The command contour is used to plot the line curve
contour(X ,Y , R4 ,[1 ,1] , "b" , "linewidth" ,2);

title( 'Stability region of  RK4')
xlabel( ' Re( h\lambda)')
ylabel( 'Im( h\lambda)')
axis square , grid on
xline(0); yline(0);





