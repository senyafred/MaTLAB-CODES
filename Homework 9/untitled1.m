% Homework 9 Problem 1
% Solving BVP  
%   y" -2y/(1+x)^2 = -4/(1+x)^2,  y(0)=0,  y(1)=1
%  by the Galerkin method.

clear all


% Set up the x-vector and the boundary conditions:
M=input('enter M = ');
h = 1/M;
x0 = 0;
xf = 1;
y0 = 0;
yf = 1;
x=[x0:h:xf];

for k=1:M-1
  x(k)= x0 + h*(k);

    for j=1:M

       %  Build the 2 subdiagonals of matrix A: a= lower; c=upper
       if k==j+1
                xmid=(x(k-1)+x(k))/2;
                Qmid= -2/(1+xmid)^2;
                c(k) = -1/h + Qmid*h/6;
       elseif k==j-1
                xmid=(x(k)+x(k+1))/2;
                Qmid= -2/(1+xmid)^2;
                a(k) = -1/h + Qmid*h/6; 
       end        
    end
end

for k=1:M
   fi= @(x)(1 + abs(x-x(k))/h);
   if k==1
       r(k) = ( 2*x(k)-4/(1+x(k))^2 )* integral(fi, x0, x(k+1));
   elseif k>1 && k<M
       r(k) = ( 2*x(k)-4/(1+x(k))^2 )*integral(fi, x(k-1), x(k+1));
   else
       r(k) = ( 2*x(k)-4/(1+x(k))^2 )*integral(fi, x(k-1), xf);
   end
    for j=1:M
        %building the main diagonal b
        if k==j
          Qmid= -2/(1+x(k))^2;
          b(k) = 2/h + Qmid*2*h/3;
        end
    end
end

const=thomas(a,b,c,r);

for j=1:M
   
    for k=1:M
        
    Y(j)= const(k)*( 1 - abs( x(k) - x(j) )/h );
    end
    
    Y(j)=Y(j)+x(j);
end

Y=[y0,Y,yf];

x=[x0,x];

for i=1:length(x)
   y(i)= 2*x(i)/( 1 + x(i)); 
end    

figure(1)
plot(x,(Y-y))
xlabel('x')
ylabel('error')
title('Error plot when h = 0.4')

figure(2)
plot(x,Y,'r--')
hold on 
plot(x,y,'g')
hold off
legend('numerical sol','exact sol')
