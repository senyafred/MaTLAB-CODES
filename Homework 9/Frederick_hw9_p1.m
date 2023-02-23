%%Problem 1

clear all;

a = 0; b=1;             %Endpoints
alpha =0;               %BC at a
beta=1;                 %BC at b

M=10;                   %Number of interval

h = (b - a)/(M+1);      %Step size

%Interval
x=[a+h:h:b-h];
x_interval=[a:h:b];

%The functions
Q = @(x) -2/(1+x)^2;
r = @(x) (2*x-4)/(1+x)^2;

for k = 1:M
  %Forming the A matrix
    for j=1:M
        A(k,j) = -j^2*pi^2*sin(j*pi*x(k)) + Q(x(k))*sin(j*pi*x(k));       
    end
    
   %The R vector
    R(k) = r(x(k));
end

%Solve the system for the c values
c = A\R'; 

%Summing the values
j = [1:M];
for k = 1:length(x)
    Z(k) = sum((c.*sin(j*pi*x(k))'));
end

%The solution
Y = Z + x;

%Updating the solution with BC
Y=[alpha,Y,beta];


% The exact solution
yexact = @(x) 2*x/(1+x);
for k = 1:length(x_interval)
    yexact1(k) = yexact(x_interval(k));
end

figure
plot(x_interval,Y, x_interval, yexact1);
xlabel('x')
ylabel('y')
legend('FEM','Exact')
title('Solution when h = 0.1')



%Finding error
error = Y - yexact1;

figure
plot(x_interval, error);
xlabel('x')
ylabel('error')
title('Error plot when h = 0.1')