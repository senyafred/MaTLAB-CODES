% The integration is done by the Modified Euler method.


% Set the parameters of the problem:
clear all;


g= 9.81;                         %Gravity
v_max=80;                        %Maximum speed(mph)
v = v_max*(1600/3600);           %Maximum speed(m/s^2)
k=g/v;                           %constant
beta= g/k;                       %constant
f = @(t,y) k*(-beta-y);          %function

t0=0;                          %Initial value of t
tf=2;                          % beginning and ending values of t
h=0.2;                         % step size in t
y0=0;                          % initial value of y

t=[t0:h:tf];                   % this is the vector of t values
%y = zeros(1,length(t));

yexact = @(t) beta*(exp(-k*t) - 1);     % Exact solution


% Solve the equation using the Modified Euler method:
y(1)=y0;
for i=1:length(t)-1
    y(i+1) = y(i) + h*f(t(i),y(i));
    y(i+1)= y(i)+h/2*(f(t(i),y(i))+f(t(i+1),y(i+1)));
end

error_ME = yexact(t)-y;   %Finding the error

plot(t,y,'b',t,yexact(t),'r','linewidth',2);  %Plotting both solutions
xlabel('Time(sec)')
ylabel('velocity')
legend(' y: Modified Euler', 'y: Exact Solution');
title('Sky diver problem', 'FontSize', 12)

figure();
plot(t,error_ME,'b--','linewidth',2);   %Plotting the error
xlabel('Time(sec)')
ylabel('Error')
title('Numerical Error', 'FontSize', 12)