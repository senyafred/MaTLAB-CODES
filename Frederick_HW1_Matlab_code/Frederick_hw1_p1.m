% This code integrates on [x0,xf] a model equation 
%   y' = a*y 
% with given constant a  and
% initial condition y(x0)=y0.
% The integration is done y the simple Euler method.


%%%%%%%%%%%%%%%
%PROBLEM 1
%%%%%%%%%%%%%%%%%

% Set the parameters of the problem:
clear all

a=1;                           %Given constant 

x0=0;                          %Initial value of x
xf=1;                          % beginning and ending values of x
h=0.1;                         % step size in x
y0=1;                          % initial value of y

x=[x0:h:xf];                   % this is the vector of x values

yexact = @ (x,y) exp(a*(x-x0))*y0;      % Exact solution 

% Solve the equation using the simple Euler method:
y(1)=y0;
for i=1:length(x)-1
    y(i+1)=y(i)+h*(a*y(i));
end

error_euler = yexact(xf)-y(end)  %Finding the error

