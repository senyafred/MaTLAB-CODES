% The integration is done by the Classical Runge Kutta method.
clear all;

g= 9.81;                         %Gravity
v_max=80;                        %Maximum speed(mph)
v = v_max*(1600/3600);           %Maximum speed(m/s^2)
k=g/v;                           %constant
beta= g/k;                       %constant

f = @(t,y) k*(-beta-y);          %differential equation

% Other parameters of the problem:
tmin = 0;
tmax = 2;                        % beginning and ending times [s]
h = 0.2;                         % step size in time

tspan = [tmin tmax];             % this is used by the solvers
t = tmin:h:tmax;

y0=0;                          %% this is the initial speed

yexact = @(t) beta*(exp(-k*t) - 1);     % Exact solution

% Solve the equation using the Modified Euler method:
[x_rk, y_rk] = Frederick_cRK(f,y0,tspan,h);

error_rk = yexact(t)-y_rk; %Finding the error

figure();
plot(t,error_rk,'b--','linewidth',2);   %Plotting the error
xlabel('Time[s]')
ylabel('Error[m/s]')
title('Numerical Error', 'FontSize', 12)