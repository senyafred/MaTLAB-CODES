clear all;

%Parameters
a = 0;  b=1;
alpha =0;       %BC at a
beta=2;         %BC at b
B1 = 1;
B2 = 2;

%h = 0.05;       %Step size
h = 0.025;
N = (b-a)/h;   %Number of iteration
x = a:h:b;     % X interval

%Setting up matrix A
upperdiag  = ones(1,N-1);
subdiag  = [ones(1,N-2),2];

%The functions
r = @(x) -4/(1+x)^2;
Q = @(x) -2/(1+x)^2;

% The exact solution
yexact = @(x) 2*x/(1+x);
for k = 1: N+1
    yexact1(k) = yexact(x(k));
end

 
for j = 2:N+1
    %% Main diagonal matrix
   if j== N+1
     a(j-1) = -2*h*(B1/B2) - (2 - h^2*Q(x(j))); 
   else
     a(j-1) = -(2 -  h^2*Q(x(j)));   
   end
    
    %% R vector
    if j == 2
         R(j-1)= h^2*r(x(j))-alpha;
    elseif j ==N+1
        R(j-1)= h^2*r(x(j))-(2*h*beta)/B2;
    else
        R(j-1)= h^2*r(x(j));
    end   
   
end
   
% Solving the system using thomas alogrithm
Y = thomas(subdiag ,a,upperdiag ,R);

%updating the solution with the left boundary condition
Y = [alpha,Y];

%Finding error
error = abs(Y - yexact1);

plot(x,error,'linewidth',1.1)
xlabel('x')
ylabel('error')
title('Error as a function of x when h = 0.025')



