%%  Problem 4

clear;

Lmax = 1;     %Maximum length
Tmax = 0.5;   %Maximum time
g_0 = 0;
g_1 = 0;
      
h= 0.1;     %Step size in x
kappa = input('enter the value for parameter k:  ');  %Step size in t

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
   for j=1:M
        if j == 1 
            u_bar(j) = g_0;
        elseif j == M
            u_bar(j) = g_1;
        else
            u_bar(j) = r*u(i,j+1) + (1-2*r)*u(i,j) + r*u(i,j-1);
        end
    end

    for j=1:M
        if j == 1             
            u(i,j) = g_0;
        elseif j == M
            u(i,j) = g_1;
        else
            u(i+1,j) = u(i,j)+(r/2)*((u(i,j+1)-2*u(i,j)+u(i,j-1)) + (u_bar(j+1) - 2*u_bar(j)+u_bar(j-1)));
        end
    end
%     figure(5000)
%     plot(x,u(i,:))
%     pause;
    
 end 
 
% % The exact solution
 yexact = @(t,x) sin(pi*x).*exp(-pi^2*t);

plot(x,u(end,:)- yexact(t(end),x))
xlabel('x')
ylabel('y')
legend('Numerical','Exact')
title('Solution when h = 0.1,  kappa = 0.006')
 

