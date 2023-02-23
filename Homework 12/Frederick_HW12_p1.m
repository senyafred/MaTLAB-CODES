%%  Problem 1

clear;

Lmax = 1;     %Maximum length
Tmax = 0.5;   %Maximum time

%Boundary condition
g_o = 0;
g_1 = 0;

h = input('enter the value for parameter h:  ');
kappa = 0.004;  %Step size in t

%Space and time intervals
t = 0:kappa:Tmax;
x = 0:h:Lmax;

M = length(x);
N = length(t);

%Initial condition
u(1,:) = sin(pi*x);

r = kappa/h^2;
 
 for i = 1:N-1
     %%Scheme
    for j = 1:M
        if j == 1 
            u(i,j) = g_o;
        elseif j == M
            u(i,j) = g_1;
        else
            u(i+1,j) = r*u(i,j+1) + (1-2*r)*u(i,j) + r*u(i,j-1);
        end
    end
    
%     figure(5000)
%     plot(x,u(i,:))
%     pause;
    
 end

%The exact solution
yexact = @(t,x) sin(pi*x).*exp(-pi^2*t);

subplot( 1, 2, 1);
plot(x,u(end,:),x, yexact(t(end),x))
xlabel('x')
ylabel('y')
legend('Numerical','Exact')
title('Solution when h = 0.05,  kappa = 0.004')

subplot( 1, 2, 2);
plot(x,u(end,:)- yexact(t(end),x))








