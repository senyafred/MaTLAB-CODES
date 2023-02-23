% This code integrates on [x0,xf] a model equation 
%   y' = sqrt(y)  
% with initial condition y(x0)=y0.
% The integration is done by the Classical Runge kutta Method.
clear all;

f = @(x,y) sqrt(y);            % function
x0=0;                          %Initial value of x
xf=1;                          % beginning and ending values of x
h=0.1;                         % step size in x
%h=0.05;
y0=1;                          % initial value of y

x=[x0:h:xf];                   % this is the vector of x values

yexact = @ (x) ((x+2)^2)/4;      % Exact solution 

% Solve the equation using the Modified Euler method:
[x_rk, y_rk] = Frederick_cRK(f,y0,x,h);


error_euler = yexact(xf)-y_rk(end); %Finding the error