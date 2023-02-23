clear all;
%Parameters
a = 2;  b=3;
alpha =-15/19;       %BC at a
beta=0;              %BC at b

h = 0.1;       %Step size
N = (b-a)/h;   %Number of interval
x = a:h:b;     % X interval

%Setting up upper and sub diagonals
upperdiag  = [2,ones(1,N-2)];
subdiag  = ones(1,N-1);

%The functions
Q = @(x) -2/x^2;
r = @(x) (x-6)/x^2;
 
for j = 1:N
    %% Main diagonal matrix
    
     a(j) = -(2 - h^2*Q(x(j)));   
    
    %% R vector
    if j == 1
         R(j)= h^2*r(x(j))+ 2*alpha*h;
    elseif j == N
        R(j)= h^2*r(x(j)) - beta;
    else
        R(j)= h^2*r(x(j));
    end   
   
end   
% Solving the system using thomas alogrithm
Y = thomas(subdiag ,a,upperdiag ,R);

%updating the solution with the boundary condition
Y = [Y,beta];

% The exact solution
yexact = @(x) (114*x-19*x^2 - 5*x^3 - 36)/(38*x);
for k = 1:length(x)
    yexact1(k) = yexact(x(k));
end

%Plotting the Solution 
figure
plot(x,Y,'k',x,yexact1,'ro','linewidth',1.7)
xlabel('x')
ylabel('Y')
legend("Numerical Solution","Exact Solution")
title('Solution when h = 0.1')

%Finding error
error = abs(Y - yexact1);

figure
plot(x,error,'linewidth',1.7)
xlabel('x')
ylabel('error')
title('Error as a function of x when h = 0.1')
