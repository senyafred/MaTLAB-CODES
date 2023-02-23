

clear all;

a = 0; b=2;
alpha =1;       %BC at a
beta=1;         %BC at b
e = 10^(-6);    %Tolerance

% Different c values for initial guess
%c = 100;             
%c = 20;         
%c = 10;           
%c = 1.5;
%c = 0.75;
%c = 0;
%c = -0.5; 
%c = -0.9;
c = -1;
%c = -1.5;

% Different c values for deep solution
%c = -7.45;             
%c = -6.95;         
%c = -6.45;           
%c = -5.95;
%c = -4;
%c = -3.9;
%c = -3.75; 
%c = -3.62;
%c = -3.6;
%c = -3.55;
%c= -3.5;
%c=-2.9;
%c=-2.5;

h = 0.1;       %Step size
%h = 0.2;       %Step size
N = (b-a)/h;   %Number of iteration

% X interval
x = a+h:h:b-h;     
x_interval=[a,x,b];

%Setting up matrix A
upperdiag  = ones(1,N-2);
subdiag  = ones(1,N-2);
%maindiag = -2*ones(1,N-1);

%Matrix A
maindiag = -2*ones(1,N-1) - c*h^2;

%The non-linear function
f = @(x,y) y.^2/(x+2);

%Intial condition
y0 = ones(1,N-1);
%y0 = piecewise_function(x);  %%Deep solution
%y0 = piecewise_function1(x);  %%Deep solution 2

%iteration count
itercount(1) = 1;

k=1;

while k
    
   if k >1000
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
    error(k)= norm(y11 - y0);
     
    %Counting the iteration
    itercount(k+1) = itercount(k) + 1;
     
      if error(k) <= e
          disp(strcat("Converges at ",string(itercount(k+1)))) 
          break
      end
      
      %Updating the initial guess
      y0 = y11;   
   
      %update k
      k = k+1;
     
end

y = [alpha,y11,beta];

%Plot solution
figure
plot(x_interval,y)
xlabel('x')
ylabel('y')
%title('solution when theta1 =-10, theta2 = -15')


%Plot Error
figure
plot(log(error))
xlabel('iteration')
ylabel('log(Error)')
title('Log(error) vs iteration (Deeper solution) when c=-6.45 and h=0.1')


function result= piecewise_function(x)
%%The piecewise function

for i = 1:length(x)
    if x(i) <= 1
        result(i) = 1 - 16*x(i);
    else
        result(i) = -31 + 16*x(i);
    end
end

end


function result= piecewise_function1(x)
%%The piecewise function 2

for i = 1:length(x)
    if x(i) <= 1
        result(i) = 1 - 20*x(i);
    else
        result(i) = -39 + 20*x(i);
    end
end

end








