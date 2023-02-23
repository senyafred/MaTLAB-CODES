%  stability regions of the P-C Method
% Z = x +1i *y 
clear all;

ZR = -2.5:0.01:2.5;   %Real
ZI = -2.5:0.01:2.5;     %Imaginary

[X,Y] = meshgrid(ZR,ZI);
z = X+1i*Y;

%Roots of quadratic equation
rho1 = 1/2*((1+z+3/4*z.^2) + sqrt((1+z+3/4*z.^2).^2 - z.^2));
rho2 = 1/2*((1+z+3/4*z.^2) - sqrt((1+z+3/4*z.^2).^2 - z.^2));

%condition
Rho1 = abs(rho1)<=1;
Rho2 = abs(rho2)<=1;

% The command contour creates a contour plot 
p1 = contour(X,Y,Rho1,[1,1], "r" , "linewidth" ,2); hold on;
p2 = contour(X,Y,Rho2,[1,1], "r" , "linewidth" ,2); hold off;

title( 'Stability region of the P-C Method')
xlabel( ' Re( h\lambda)')
ylabel( 'Im( h\lambda)')