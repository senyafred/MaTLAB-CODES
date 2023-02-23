%%Solving x^2u" -2u = x - 6     u'(2) = -15/19 u(3) = 0 
%using the Shooting method

clear all;

%%Boundary condition
alpha = -15/19;  beta=0;

x0=2;                          %Initial value of x
xf=3;                          % final Value of x
h=0.1;                         % step size 

x=[x0 xf];                     % this is the vector of x values

% yint=[1 0];                    %Initial conditions for 1st Auxillary
% yint2=[0 alpha];               %Initial conditions for 2nd Auxillary

yint=[0 alpha];                    %Initial conditions for 1st Auxillary
yint2=[1 0];               %Initial conditions for 2nd Auxillary

% Solving Auxilliary IVPs using the Modified Euler method:
[x1,u1]=Frederick_Test2_modified_Euler(@Auxi1,x,h,yint);

[x2,v2]=Frederick_Test2_modified_Euler(@Auxi2,x,h,yint2);
 
 %Finding parameter theta
 theta = (beta - u1(1,end))/v2(1,end);
 
 %The BVP solution
 y= u1 + theta*v2;     
 
 % The exact solution
yexact = @(x) (114*x-19*x^2 - 5*x^3 - 36)/(38*x);
for k = 1:length(x1)
    yexact1(k) = yexact(x1(k));
end
 
 %Plotting the Solution 
plot(x1,y(1,:),'b',x1,yexact1,'ro','linewidth',1.7)
xlabel('x')
ylabel('Y')
legend("Numerical Solution","Exact Solution")
title('Solution using shooting method when h = 0.1')

%Finding error
error = abs(y(1,:) - yexact1);

figure
plot(x1,error,'linewidth',1.7)
xlabel('x')
ylabel('error')
title('Error as a function of x when h = 0.1')

function dy = Auxi1(x, u)
%First Auxilliary Equation
  dy = [u(2); (2/x^2)*u(1)+((x-6)/x^2)];
end

function dy = Auxi2(x, v)
%Second Auxilliary Equation
  dy = [v(2); (2/x^2)*v(1)];
end