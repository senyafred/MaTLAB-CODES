%  stability regions of the 2nd-order Adams-Bashforth Method
% Z = x +1i *y 
clear all;

ZR = -1.5:0.005:1.5;   %Real
ZI = -1.5:0.005:1.5;     %Imaginary

[X,Y] = meshgrid(ZR,ZI);
z = X+1i*Y;

%Roots of quadratic equation
rho1 = 1/2*((1+3/2*z) + sqrt((1+3/2*z).^2 - 2*z));
rho2 = 1/2*((1+3/2*z) - sqrt((1+3/2*z).^2 - 2*z));

%condition
Rho1 = abs(rho1)<=1;
Rho2 = abs(rho2)<=1;

% The command contour creates a contour plot 
p1 = contour(X,Y,Rho1,[1,1], "k" , "linewidth" ,2); hold on;
p2 = contour(X,Y,Rho2,[1,1], "k" , "linewidth" ,2); hold off;

title( 'Stability region of 2nd-order Adams-Bashforth Method')
xlabel( ' Re( h\lambda)')
ylabel( 'Im( h\lambda)')
