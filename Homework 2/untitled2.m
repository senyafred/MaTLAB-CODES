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
%
% Examples:
 [y, t]=rkf45(@frederick_fun4_hw2_p3,0,4,0,0.2,1e-3);         % here 'myfunc' is a user-defined function in M-file
 
 
 % % The exact solution
yexact1 = @(t) (g/k1)*(exp(-k1*t) - 1); 
yexact2 = @(t) ((g + y02*k2)*exp(-k2*(t-t02)) - g)/k2;

x1 = 0:h:2;
x2 = 2.2:h:4;

exact=[yexact1(x1),yexact2(x2)];
 
figure(10501);
plot(t,y,'k','linewidth',2)%Plotting both solutions
xlabel('Time[s]')
ylabel('velocity[m/s]')
legend(' y: cRK', 'y: Exact Solution');
title('Sky diver problem', 'FontSize', 12) 


 
 
% [w t]=rkf45(inline('sin(y*t)','t','y'),0,1,1,0.1,0.01);
% f=inline('sin(y(1))-cos(y(2))','t','y');
% [y t]=rkf45(f,0,1,1,0.01,1e-5);

function [y, tt] = rkf45(f,a,b,ya,h,rtol)

% Compute the constants once
c30 = 3/8;
c31 = 3/32;
c32 = 9/32;
c40 = 12/13;
c41 = 1932/2197;
c42 = -7200/2197;
c43 = 7296/2197;
c51 = 439/216;
c52 = -8;
c53 = 3680/513;
c54 = -845/4104;
c61 = -8/27;
c62 = 2;
c63 = -3544/2565;
c64 = 1859/4104;
c65 = -11/40;
cw1 = 25/216;
cw3 = 1408/2565;
cw4 = 2197/4104;
cw5 = -1/5;
cz1 = 16/135;
cz3 = 6656/12825;
cz4 = 28561/56430;
cz5 = -9/50;
cz6 = 2/55;
ce1 = 1/360;
ce3 = -128/4275;
ce4 = -2197/75240;
ce5 = 1/50;
ce6 = 2/55;

% Absolute tolerance
atol = 1e-13;
alpha = 0.8;
k = 0;
% Initial time moment
i = 1;
tt(1) = a;
t = a;
% Initial condition
y(1,:) = ya;
wi = ya;
% If it is the last iteration, then lastit = 1, otherwise lastit = 0
lastit = 0;
while lastit == 0
    % Stretch the step if within 10% of b-t
    if t + 1.1*h > b
        h = b - t;
        lastit = 1;
    end
   
    % Compute the step
    s1 = f(t,wi);
    s2 = f(t + 0.25 * h, wi + 0.25*h*s1);
    s3 = f(t + c30 * h, wi + c31 * h * s1 + c32 * h * s2);
    s4 = f(t + c40 * h, wi + c41 * h * s1 + c42 * h * s2 + c43 * h * s3);
    s5 = f(t + h, wi + c51 * h * s1 + c52 * h * s2 + c53 * h * s3 + c54 * h * s4);
    s6 = f(t + 0.5 * h, wi + c61 * h * s1 + c62 * h * s2 + c63 * h * s3 + c64 * h * s4 + c65 * h * s5);
    w = wi + h * (cw1 * s1 + cw3 * s3 + cw4 * s4 + cw5 * s5);
    z = wi + h * (cz1 * s1 + cz3 * s3 + cz4 * s4 + cz5 * s5 + cz6 * s6);
    e = h * norm(ce1 * s1 + ce3 * s3 + ce4 * s4 + ce5 * s5 + ce6 * s6);
    
    % Target tolerance for this step
    T = rtol * norm(wi) + atol;
    if e <= T % In case the tolerance is met
        t = t + h;
        h = alpha*h*(T/e)^0.2;
        i = i + 1;
        tt(i) = t;
        wi = z;
        y(i,:) = z;
        k = 0;
    elseif k == 0 % Tolerance is not met for the first time in this step
        h = alpha*h*(T/e)^0.2;
        k = k + 1;
        lastit = 0;
    else % Tolerance is not met more than once in this step
        h = h / 2;
        lastit = 0;
    end
        
end
end