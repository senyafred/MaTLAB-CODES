% The integration is done by ODE-Solver 0de45.

clear all;

global g k1 k2

% parameters of the problem:
vmax_mph1 = 80;                    % max speed in miles per hour
vmax_mph2 = 4;
vmax1 = vmax_mph1*(1600/3600);     % max speed in meters per second
vmax2 = vmax_mph2*(1600/3600);
g = 9.8;                           % acceleration of gravity in m/s^2
k1 = g/vmax1;                      % friction constants
k2 = g/vmax2;


% Other parameters of the problem:
tmin = 0;                        % beginning and ending times [s]
tmax = 4;                        
h = 0.2;                         % step size in time

t = tmin:h:tmax;                 % this is used by the solvers
t02 = 2;                         % initial time for second ODE
y02 = -15.1;                     % initial speed for second ODE


y0=0;                          %% this is the initial speed

% Solve the equation using ode45:
[t_ode,y_ode] = ode45(@frederick_fun4_hw2_p3,t,y0); 

% The exact solution
yexact1 = @(t) (g/k1)*(exp(-k1*t) - 1);  
yexact2 = @(t) ((g + y02*k2)*exp(-k2*(t-t02)) - g)/k2;

x1=0:h:2;
x2= 2.2:h:4;

exact=[yexact1(x1),yexact2(x2)];    %Patching the solution

% Find the error
error_ode = y_ode' - exact;

% Save the error to file:
save('frederick_hw2_p4','error_ode')

figure(10501);
plot(t,exact,'b--',t_ode,y_ode,'k','linewidth',2)%Plotting both solutions
xlabel('Time[s]')
ylabel('velocity[m/s]')
legend(' y: ODE45', 'y: Exact Solution');
title('Sky diver problem', 'FontSize', 12)
