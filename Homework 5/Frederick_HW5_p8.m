% This code integrates on [x0,xf] a model equation
% Set the parameters of the problem:
clear all;

h = input('Enter the step size: h =');  
alpha = input('Enter the friction coefficent: alpha =');
xf = input('Enter maximum time: t =');

f1 = @(x,y2) y2;              % system of ODE where y2 =v
f2 = @(x,y1,y2,w_o,alpha) -w_o^2*y1 -2*alpha*y2; 

x0=0;                          %Initial value of x
                                 
w_o = 1;                        %damped oscillator

x=x0:h:xf;                   % this is the vector of x values

w = sqrt(w_o^2 - alpha^2);

y1(1)=0;   %initial conditions 
y2(1)=1;

z1(1)=0;    %initial conditions 
z2(1)=1;
%Solve the equation using the Regular Euler
for i=1:length(x)-1
    y1(i+1) = y1(i) + h*f1(x(i),y2(i));
    y2(i+1) = y2(i) + h*f2(x(i),y1(i),y2(i),w_o,alpha);
    
end
E1 = (y2.^2 + w^2*y1.^2)/2;  

H1 = 0.112414 -E1;  %Hamiltonian Error

%Solve the equation using the Sympletic Euler
for i=1:length(x)-1
    z1(i+1)= z1(i)+h*(f1(x(i),z2(i)));
    z2(i+1)= z2(i)+h*(f2(x(i),z1(i+1),z2(i),w_o,alpha));
    
end
E2 = (z2.^2 + w^2*z1.^2)/2;

H2 = 0.112414 -E2;  %Hamiltonian Error

% The exact solution
yexact = @(t) (1/w)*(exp(-alpha*t).*sin(w*t));  
yexact1 = @(t) (exp(-alpha*t).*(w.*cos(w*t) - alpha.*sin(w*t))).*(1/w); 

%%Subploting the graph
subplot(1,3,1);
plot(yexact(x),yexact1(x),y1,y2)
subplot(1,3,2);
plot(yexact(x),yexact1(x),z1,z2)
subplot(1,3,3)
plot(x,H1); hold on
plot(x,H2);






