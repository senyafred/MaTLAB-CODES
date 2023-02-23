
clear all;

a = 0;
b=2;
alpha =1;       %BC at a
beta=1;         %BC at b
e = 10^(-6);    %Tolerance

%h = 0.1;       %Step size
h = 0.2;       %Step size
N = (b-a)/h;   %Number of iteration
x = a:h:b;     % X interval

%Setting up matrix A
upperdiag  = ones(1,N-2);
subdiag  = ones(1,N-2);

%The non-linear function
f = @(x,y) y.^2/(x+2);
df = @(x,y) 2*y/(x+2);

%Intial condition
%y0 = ones(1,N+1);          %%Shallow IC
y0 = piecewise_function(x);  %%Deep IC

%Error conditions
e0 = 0;  eN=0;

%itercount(1) = 1;
k = 1;

while k
      
  if k >1000
       disp("iteration do NOT converge!!!!!!")
       break
  end
    
  %Main diagonal
    for j = 2:N
        a(j-1) = -(2 + h^2*df(x(j),y0(j)));
    end
    
    %R vector
    for j = 2:N
        if j ==2
            r(j-1)= y0(j+1) -2*y0(j) + y0(j-1) - h^2*f(x(j),y0(j))-e0;
        elseif j == N
            r(j-1)= y0(j+1) -2*y0(j) + y0(j-1) - h^2*f(x(j),y0(j))-eN;
        else
            r(j-1)= y0(j+1) -2*y0(j) + y0(j-1) - h^2*f(x(j),y0(j));
        end           
    end
     
    %Solving the system using Thomas algorithm
    discrep = thomas(subdiag ,a,upperdiag ,r);
    
    %New solution
    y1 = y0(2:N) - discrep;
    
    %Update solution with BC
    y = [alpha,y1,beta];
    
    %  Find the error
    error(k) = norm(y - y0);
    
       if error(k) <= e   
           disp(strcat("Converges at ",string(k)))
           break
       end
       
      y0 = y;
      
     %update the iteration
     k = k+1;
      
end

figure
plot(x,y)

figure
plot(log(error)); 
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

