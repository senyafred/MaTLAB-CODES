%%Problem 7

clear all;

a = 0; b=2;     %endpoints
alpha =1;       %BC at a
beta=1;         %BC at b
e = 10^(-6);    %Tolerance

%h = 0.1;       %Step size
h = 0.2;
N = (b-a)/h;   %Number of iteration

% X interval
x = a+h:h:b-h;     
x_interval=[a,x,b];

%Setting up matrix A
upperdiag  = ones(1,N-2);
subdiag  = ones(1,N-2);
maindiag = -2*ones(1,N-1);

%The non-linear function
f = @(x,y) y.^2/(x+2);

%Intial condition
y0 = ones(1,N-1);     %%Shallow condition
%y0 = piecewise_function(x);  %%Deep solution

%iteration count
k = 1;

while k
    
  if k >1000
       disp("iteration do NOT converge!!!!!!")
       break
  end
    
   % the R vector
    for j = 1:N-1
        if j ==1
            r(j)= h^2*f(x(j),y0(j))-alpha;
        elseif j == N-1
            r(j)= h^2*f(x(j),y0(j))-beta;
        else
            r(j)= h^2*f(x(j),y0(j));
        end           
    end
    
    %Solve the system using Thomas algorithm
     y11 = thomas(subdiag ,maindiag, upperdiag ,r);
     
     %FInd the error
     error(k)= norm(y11 - y0);
      
      if error(k) <= e
          disp(strcat("Converges at ",string(k)))
          break
      end
      
     y0 = y11;
     
     %update the iteration
     k = k+1;
end

%update solution with boundary conditions
y = [alpha,y11,beta];

%Plot solution
plot(x_interval,y)
xlabel('x')
ylabel('y')
title('y as a function of x (Deep Solution) when h = 0.2')

%Plot the error
figure
plot(log(error),'linewidth',1.1)
xlabel('iteration')
ylabel('log(error)')
title('Maximum discrepancy (Deep Solution) when h = 0.2')

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


