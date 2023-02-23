% Test 2 Problem 4
% Solving BVP  
%  u_xx + (6sech^2(x) - lambda)u - u^3 = 0
%  by the Modified Picard method.

clear all;

a = -15; b=15;
alpha =0;       %BC at a
beta=0;         %BC at b
tol = 10^(-9);    %Tolerance

% c values for initial guess
c = input('enter the value for parameter C:  C = ');
lamda = input('enter the value for parameter {\lambda}:  {\lambda} = ');

h = 0.05;       %Step size
N = (b-a)/h;    %Number of iteration
max_Niterations =1000;   %Maximum number of iteration

% X interval
x = a+h:h:b-h;     
x_interval=[a,x,b];

%Setting up upper and lower vectors
upperdiag  = ones(1,N-2);
subdiag  = ones(1,N-2);

%Main diagonal
maindiag = -2*ones(1,N-1) - c*h^2;

%The non-linear function
f = @(x,y) -(6*sech(x)^2 - lamda)*y + y^3;

%Intial conditions
y1 = @(x) exp(-x^2);      %For lambda = 2.25

%y1 = @(x) x*exp(-x^2);     %For lambda = 0.778
for k = 1:length(x)
    y0(k) = y1(x(k));
end

%iteration count
k=1;

%Assign initial error
error(1) = 1;

% # of iterations after which the solution is plotted within the loop below
Nplot = 20;

while error(k) > tol && k < max_Niterations
    
   if k >max_Niterations
       disp("iteration do NOT converge!!!!!!")
       break
   end
       
    for j = 1:N-1
        if j ==1
            r(j)= h^2*f(x(j),y0(j))-alpha;
        elseif j == N-1
            r(j)= h^2*f(x(j),y0(j))-beta;
        else 
            r(j)= h^2*f(x(j),y0(j));
        end
    end
    
    %R vector 
    r_ = r - c*h^2*y0;
    
    %using thomas algorithm
    y11 = thomas(subdiag ,maindiag, upperdiag ,r_);
    
    %Checking for error
    error(k+1)= norm(y11 - y0,inf);
     
      if error(k) <= tol 
          break
      end
      
      %Updating the initial guess
      y0 = y11;   
   
     % Plot solutions after Nplot iterations:
    if rem(k, Nplot) == 0
        figure(80902);
        plot(x_interval, [alpha y0 beta]);
        xlabel('X','fontsize',16);  
        ylabel('Y','fontsize',16); 
        title(['iteration # ' int2str(k)], 'fontsize', 16)
        disp('Hit any key to continue with iterations.')
        pause
    end
     
      %Counting the iteration
      k = k+1;
end

%Number of iterations
sprintf('converges in %d  iterations', k)

%Update the solution with BC
y = [alpha,y11,beta];

% %Plot solution
subplot( 1, 2, 1);
plot(x_interval,y)
xlabel('x')
ylabel('y')

%Plot Error
subplot( 1, 2, 2);
plot([1:k], log10(error),'r');
%plot(log10(error))
xlabel('iteration`s number','fontsize',16);  
ylabel('log_{10} of error','fontsize',16); 
title('Error when c= 2.6 and lambda = 2.25')

