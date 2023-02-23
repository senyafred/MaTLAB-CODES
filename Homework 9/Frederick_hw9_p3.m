%%Problem 3

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
r = @(x) (2*x-4)/(1+x).^2;

for k = 1:M-1
    %Upper and Lower diagonals
    for j=1:M     
        if k == j+1
            upperdiag(k) =  1/h + (h/6)*Q(x(j)+h/2);
        elseif k == j-1 
             subdiag(k) =  1/h + (h/6)*Q(x(k)+h/2);
        else
            continue
        end
    end   
end

for k = 1:M
    %Main diagonal 
    for j=1:M     
        if k == j
            maindiag(k) =  -2/h + (2*h/3)*Q(x(j));
        end
    end
end


for k=1:M
    %The R function
   if k==1
       R(k) = integral(@(x) r(x).*(1 + abs(x-x(k))/h), a , x(k+1));
       
   elseif k==M
       R(k) = integral(@(x) r(x).*(1 + abs(x-x(k))/h), x(k-1) ,b );
       
   else
       R(k) = integral(@(x) r(x).*(1 + abs(x-x(k))/h), x(k-1), x(k+1));
   end
end

%Solve the system for the c values
c = thomas(subdiag,maindiag,upperdiag,R);

%Summing the values
j = [1:M];
for k = 1:length(x)
    Z(k) = sum((c.*(1-abs(x(j)-x(k))/h)));
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