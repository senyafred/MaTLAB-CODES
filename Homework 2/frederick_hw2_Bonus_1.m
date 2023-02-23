clear all;

global g k1 k2

% parameters of the problem:
vmax_mph1 = 80;                   % max speed in miles per hour
vmax_mph2 = 4;
vmax1 = vmax_mph1*(1600/3600);     % max speed in meters per second
vmax2 = vmax_mph2*(1600/3600);
g = 9.8;                         % acceleration of gravity in m/s^2
k1 = g/vmax1;                     % friction constant
k2 = g/vmax2;


% Other parameters of the problem:
tmin = 0;
tmax = 4;                        % beginning and ending times [s]
h = 0.2;                         % step size in time

tspan = [tmin tmax];             % this is used by the solvers
%t = tmin:h:tmax;
t02 = 2;                         % initial time for second ODE
y02 = -15.1;                     % initial speed for second ODE

y0=0;                          %% this is the initial speed

% Solve the equation using the Classical Rk method:
[t_rkf, Y_rkf] = Frederick_RKF(@frederick_fun4_hw2_p3);


