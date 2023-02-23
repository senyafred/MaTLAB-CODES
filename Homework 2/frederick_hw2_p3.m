% The integration is done by the Classical Runge Kutta method.

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

tspan = [tmin tmax];             % this is used by the solvers
t02 = 2;                         % initial time for second ODE
y02 = -15.1;                     % initial speed for second ODE

y0=0;                            %% this is the initial speed

% Solve the equation using the Classical Rk method:
[t_rk, Y_rk] = Frederick_cRK(@frederick_fun4_hw2_p3,y0,tspan,h);
disp(Y_rk(end))

% The exact solution
yexact1 = @(t) (g/k1)*(exp(-k1*t) - 1); 
yexact2 = @(t) ((g + y02*k2)*exp(-k2*(t-t02)) - g)/k2;

x1 = 0:h:2;
x2 = 2.2:h:4;

exact=[yexact1(x1),yexact2(x2)];   %%Patching the two solutions

error_rk = Y_rk-exact;             %Finding the error

% Save the error to file:
save('frederick_hw2_p3','error_rk')

figure(10501);
plot(t_rk,exact,'b--',t_rk,Y_rk,'k','linewidth',2)%Plotting both solutions
xlabel('Time[s]')
ylabel('velocity[m/s]')
legend(' y: cRK', 'y: Exact Solution');
title('Sky diver problem', 'FontSize', 12)

figure(10502);
plot(t_rk,error_rk,'r--','linewidth',2)   %Plotting the error   
xlabel('Time[s]')
ylabel('Error[m/s]')
title('Numerical Error', 'FontSize', 12)
