% This code integrates on [x0,xf] a model equation 
%   y' = -20y  
% with initial condition y(x0)=y0.

clear all;

k0 = -20;                       %Constant from differential equation
f = @(x,y) k0*y;                % function

x0=0;                            %Initial value of x
xf=1.5;                          % beginning and ending values of x
h=0.125;                         % step size in x


y0=1;                          % initial value of y

x=[x0:h:xf];                   % this is the vector of x values

yexact = @(x) y0*exp(k0*(x-x0));      % Exact solution 

% Solve the equation using the Simple Euler method:
[x_s, y_s] = Frederick_HW4_p7_Euler(f,y0,x,h);

% Solve the equation using the classical Runge Kutta method:
[x_rk, y_rk] = Frederick_HW4_p7_cRK(f,y0,x,h);

% Solve the equation using the Implicit Euler method:
[x_I, y_I] = Frederick_HW4_p7_Implicit(y0,x,h);


figure(10501);
plot(x_s,y_s,'b','linewidth',2);  %Plotting both solutions
xlabel('x')
ylabel('y')
title('Simple Euler Method: dy/dx= -20y', 'FontSize', 12)

figure(10502);
plot(x_rk,y_rk,'k',x_I,y_I,'b--',x,yexact(x),'r+-','linewidth',2);  %Plotting both solutions
xlabel('x')
ylabel('y')
legend(' cRK', 'Implicit Euler','Exact Solution');
title('cRK Method Vs Implicit Euler Method Vs Exact Solution: dy/dx= -20y', 'FontSize', 12)


