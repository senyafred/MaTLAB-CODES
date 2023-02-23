% This code integrates on [x0,xf] the simple harmonic oscillator model
%   

clear all;

x0=0;                            %Initial value of x
xf=20;                          % ending value of x
h=0.1;                          % step size in x


y0 = 0;
y1 = 1;                         % initial value chosing for y2

x=[x0:h:xf];                   % this is the vector of x values

% Solve the equation using the Implicit Euler method:
[x_I, y1_I,y2_I] = Frederick_Test1_p4_Implicit(y0,y1,x,h);


figure(10502);
plot(y1_I,y2_I,'r');  %Phase plot
xlabel('y')
ylabel('v')
title('Phase Portrait of the simple oscillator model')