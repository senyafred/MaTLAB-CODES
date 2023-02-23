% This code integrates on [x0,xf] a model equation 
%   y' = sqrt(y)
% with initial condition y(x0)=y0.
% The integration is done by the Midpoint method.


% Set the parameters of the problem:
clear all;

f = @(x,y) sqrt(y);            % function
x0=0;                          %Initial value of x
xf=1;                          % beginning and ending values of x
%h=0.1;                         % step size in x
h=0.05;
x=[x0:h:xf];                   % this is the vector of x values
y0=1;                          % initial value of y

% Solve the equation using the midpoint method:

halfh = h / 2;
y(1) = y0;
for i = 1 : length(x)-1
    z = y(i) + halfh * f(x(i),y(i));
    y(i+1) = y(i) + h * f(x(i)+halfh,z);
end

yexact = @ (x) ((x+2)^2)/4;      % Exact solution 

error_euler = yexact(xf)-y(end)  %Finding the error

