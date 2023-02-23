%%Plotting the global errors

clear all;

t = [0:0.2:4];

load('frederick_hw2_p3');   %Load error from cRK

load('frederick_hw2_p4');   %Load error from ODE45

figure(10501);
plot(t, error_rk,'r-',t,error_ode,'b--','linewidth',2)  %Plotting both solutions
xlabel('Time[s]')
ylabel('Error[m/s]')
legend(' Error: cRK', 'Error: Ode45');
title('Errors of the Numerical Solutions in Problem 3 & 4', 'FontSize', 12)