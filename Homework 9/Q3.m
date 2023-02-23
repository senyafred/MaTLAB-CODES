
%%Problem 7

clear all;

a = 0; b=1;     %endpoints
alpha =0;       %BC at a
beta=1;         %BC at b

M=10;



% M1 = length(x);
% for i=1:M1-1,
%   h1(i) = x(i+1)-x(i);
% end
h = (b - a)/(M+1);       %Step size
x=[a:h:b];




%The function
Q = @(x) -2/(1+x)^2;

r = @(x) (2*x-4)./(1+x).^2;



for k = 2:M+1
    for j=1:M     
        if k-1 == j
            A(k-1,j) =  -2/h + (2*h/3)*Q(x(j));
        elseif k-1 == j+1
            A(k-1,j) =  1/h + (h/6)*Q(x(j)+h/2);
        elseif k-1 == j-1 
             A(k-1,j) =  1/h + (h/6)*Q(x(k)+h/2);
        else
            A(k-1,j) = 0;
        end
    end
    
  
 R(k-1) = integral(@(x) r(x).*(1 + (x-x(k))/h), x(k-1), x(k)) + integral(@(x) r(x).*(1 - (x-x(k))/h), x(k), x(k+1));
    
    
end


c = A\R'; 

for k = 1:length(x) 

    for j = 1: M
    if j ==1
        Z(k) = c(j).*(1 - abs(x(k)-x(j+1))./h);
    else
        Z(k) = Z(k)+ c(j).*(1 - abs(x(k)-x(j+1))./h);
    end 
    end
    
end

Z = Z';

Y = Z + x';

% The exact solution
yexact = @(x) 2*x/(1+x);
for k = 1:length(x)
    yexact1(k) = yexact(x(k));
end

figure
plot(x,Y, x, yexact1);

%Finding error
error = abs(Y - yexact1);

figure
plot(x, error);
xlabel('x')
ylabel('error')
title('Error plot when h = 0.1')



